library(MASS)
library(ggplot2)
library(rethinking)

# Simulate data for the baseline model
simulated.baseline <- function(sp=25, sites=300, Dbeta=NULL, Dgamma=NULL,
                               beta_s = 0.3, beta_mu = 0, gamma_s=0.1, gamma_mu=0,
                               beta_nu = 2, beta_rho=4, gamma_nu=1, gamma_rho=10,
                               alpha_mu = -0.5, alpha_s=1
                               ){

        # Environmental predictor
        X <- rnorm(sites)
        
        # Generate basic distance matrix D1
        if(is.null(Dbeta)){
            Dbeta <- (as.matrix(dist(1:sp))/sp)
        }
        # Generate basic distance matrix D2
        if(is.null(Dgamma)){
            Dgamma <- (as.matrix(dist(1:sp))/sp)
        }
        
        # Generate variance-covariance structures
        Sigma_beta <- beta_nu*exp(-beta_rho*(Dbeta**2)) + diag(sp)*beta_s
        Sigma_gamma <- gamma_nu*exp(-gamma_rho*(Dgamma**2)) + diag(sp)*gamma_s
        
        # Sample beta and gamma
        beta <- mvrnorm(mu = rep(beta_mu, times = sp), Sigma = Sigma_beta)
        gamma <- exp(mvrnorm(mu = rep(gamma_mu, times = sp), Sigma = Sigma_gamma))

        # Sample alphas
        alpha <- exp(rnorm(sp, alpha_mu,alpha_s))

        # Generate dataset object
        dataset <- expand.grid(site=1:sites, id=1:sp)
        
        # Populate the object
        dataset$S1 <- X[dataset$site]
        dataset$beta <- beta[dataset$id] 
        dataset$alpha <- alpha[dataset$id]
        dataset$gamma <- gamma[dataset$id]
        
        # Calculate probability
        dataset$p <- exp(-dataset$alpha - dataset$gamma*(dataset$beta - dataset$S1)**2)
        
        # Sample occurrences
        dataset$obs <- rbinom(n = length(dataset$S1), size = 1, prob = dataset$p)
        
        # check if all species have at least a couple of occurrences (in the article we set this minimum number to be 20)
        test <- colSums(matrix(dataset$obs, sites, sp))
        if(any(test<2)){
                warning('There is at least one species that have less than two occurrences. You might want to drop those sp or rerun (try different parameter values)')
        }
        
        # Return object with all the data, parameter values and distance matrices
        return(list(dataset=dataset, Dbeta=Dbeta, Dgamma=Dgamma))
}

# Simulate data for the baseline model
simulated.generror <- function(sp=25, sites=300, Dbeta=NULL, Dgamma=NULL,
                               beta_s = 0.3, beta_mu = 0, gamma_s=0.1, gamma_mu=0,
                               beta_nu = 2, beta_rho=4, gamma_nu=1, gamma_rho=10,
                               alpha_mu = -0.5, alpha_s=1,
                               isnu=T, islambda=T, 
                               nu_mu=0.5, nu_s=0.1, lambda_mu=1.1, lambda_s=0.4
){
        
        # Environmental predictor
        X <- rnorm(sites)
        
        # Generate basic distance matrix D1
        if(is.null(Dbeta)){
                Dbeta <- (as.matrix(dist(1:sp))/sp)
        }
        # Generate basic distance matrix D2
        if(is.null(Dgamma)){
                Dgamma <- (as.matrix(dist(1:sp))/sp)
        }
        
        # Generate variance-covariance structures
        Sigma_beta <- beta_nu*exp(-beta_rho*(Dbeta**2)) + diag(sp)*beta_s
        Sigma_gamma <- gamma_nu*exp(-gamma_rho*(Dgamma**2)) + diag(sp)*gamma_s
        
        # Sample beta and gamma
        beta <- mvrnorm(mu = rep(beta_mu, times = sp), Sigma = Sigma_beta)
        gamma <- exp(mvrnorm(mu = rep(gamma_mu, times = sp), Sigma = Sigma_gamma))
        
        # Sample alphas, nu, and lambda
        alpha <- exp(rnorm(sp, alpha_mu,alpha_s))
        
        if(isnu){
                # Parameter controlling the tails
                nu <- abs(rnorm(sp,nu_mu, nu_s))+1
        }else{
                # Gaussian tails
                nu <- rep(2,sp) 
        }
        
        if(islambda){
                # Parameter controlling the skew
                lambda <- inv_logit(rnorm(sp, lambda_mu,lambda_s))*2-1
        }else{
                # No skew
                lambda <- rep(0,sp)
        }
        
        # Transform gamma and beta
        gamma_v <- gamma * sqrt((pi*(1+3*lambda*lambda)*base::gamma(3/nu)-(16**(1/nu))*lambda*lambda*base::gamma(0.5+1/nu)*base::gamma(0.5+1/nu)*base::gamma(1/nu))/(pi*base::gamma(1/nu)))
        beta_m <- beta - 2**(2.0/nu)*lambda*base::gamma(0.5+1.0/nu)/(sqrt(pi)*gamma_v)

        # Generate dataset object
        dataset <- expand.grid(site=1:sites, id=1:sp)
        
        # Populate the object
        dataset$S1 <- X[dataset$site]
        dataset$beta <- beta[dataset$id] 
        dataset$alpha <- alpha[dataset$id]
        dataset$gamma <- gamma[dataset$id]
        dataset$lambda <- lambda[dataset$id]
        dataset$nu <- nu[dataset$id]
        dataset$betam <- beta_m[dataset$id]
        dataset$gammav <- gamma_v[dataset$id]
        
        # Calculate probability
        dataset$p <- exp(-dataset$alpha - (dataset$gammav*abs(dataset$S1-dataset$betam)/(1+dataset$lambda*sign(dataset$S1-dataset$betam)))**dataset$nu)
        
        # Sample occurrences
        dataset$obs <- rbinom(n = length(dataset$S1), size = 1, prob = dataset$p)
        
        # check if all species have at least a couple of occurrences (in the article we set this minimum number to be 20)
        test <- colSums(matrix(dataset$obs, sites, sp))
        if(any(test<2)){
                warning('There is at least one species that have less than two occurrences. You might want to drop those sp or rerun (try different parameter values)')
        }
        
        # Return object with all the data, parameter values and distance matrices
        return(list(dataset=dataset, Dbeta=Dbeta, Dgamma=Dgamma))
}

# Visualize simulated data
plot.distributions <- function(dataset){
        ggplot(data=dataset, aes(x=S1, y=p, colour=as.factor(id)))+
                geom_line()+
                scale_x_continuous(expand = expansion(add = c(0, 0)))+
                theme_bw()+
                theme(legend.position = "none",
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank())
}

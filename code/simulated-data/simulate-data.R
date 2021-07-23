library(MASS)
library(ggplot2)

# Simulate data for the baseline model
simulated.baseline <- function(sp=25, sites=300, Dbeta=NULL, Dgamma=NULL,
                               beta_s = 0.3, beta_mu = 0, gamma_s=0.1, gamma_mu=0,
                               beta_nu = 3, beta_rho=4, gamma_nu=1, gamma_rho=10,
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
                warning('There is at least one species that have less than two occurrences. You might want to rerun (try different parameter values)')
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



# Generate fake data to test the extent to which the model works
simulated.generr.data <- function(){
        
        set.seed(30)
        # Define system dimensions
        N <- 20
        sites <- 300
        
        # Environmental predictors for each site
        e1 <- rnorm(sites)
        e2 <- rnorm(sites)
        
        # uncorrelated coefficients for each species and parameters
        sigma1 <- 0.3 # sd beta1
        mean1 <- -1 # mean beta1
        sigma2 <- 0.4 # sd beta2
        mean2 <- 1.5 # mean beta2
        rho <- 0.4 # correlation between betas
        
        # Coefficients for generating the variance-covariance matrix
        nu <- 3
        s <- 0.5
        Dis <- (as.matrix(dist(1:N))/N)
        
        # coefficients for each species
        Sigma <- nu*exp(-1/(s*s)*(Dis^2)) + diag(N)*sigma1
        z1 <- mvrnorm(mu = rep(0, times = N), Sigma = Sigma)
        
        # Generate correlations
        beta1 <- z1
        
        vec <- c(1:round(N*0.5), 1:round(N*0.5))
        Dis_sigma <- as.matrix(dist(vec))+1-diag(N)
        Dis_sigma <- (Dis_sigma/max(Dis_sigma))
        
        Sigma <- 1*exp(-1/(0.3*0.3)*(Dis_sigma^2)) + diag(N)*0.1
        sigma_beta1 <- exp(mvrnorm(mu = rep(-1, times = N), Sigma = Sigma))
        # Sigma <- 1*exp(-1/(0.2*0.2)*(Dis^2)) + diag(N)*0.1
        # alpha <- exp(mvrnorm(mu = rep(0, times = N), Sigma = Sigma))
        alpha <- exp(rnorm(N, 0,1))
        nu <- exp(rnorm(N,-1,0.5))+1
        sigma_beta1 <- sqrt( (sigma_beta1 * gamma(3/nu)) / gamma(1/nu))
        
        # Simulate data
        dataset <- expand.grid(site=1:sites, id=1:N)
        dataset$S1 <- e1[dataset$site]
        dataset$S2 <- e2[dataset$site]
        dataset$beta1 <- beta1[dataset$id] 
        dataset$alpha <- alpha[dataset$id]
        dataset$nu <- nu[dataset$id]
        dataset$sigma_beta1 <- sigma_beta1[dataset$id]
        
        dataset$p <- exp(-alpha[dataset$id] - sigma_beta1[dataset$id]*abs(beta1[dataset$id] - dataset$S1)**nu[dataset$id])
        
        dataset$obs <- rbinom(n = length(dataset$S1), size = 1, prob = dataset$p)
        dataset <- data.frame(id=dataset$id, p=dataset$p, obs=dataset$obs, alpha=dataset$alpha, beta1=dataset$beta1, sigma_beta1=dataset$sigma_beta1,nu=dataset$nu, S1=dataset$S1, S2=dataset$S2)                
        return(list(dataset=dataset, corr=Dis, corr2=Dis_sigma, corr3=Dis))
}

# Generate fake data to test the extent to which the model works
simulated.skew.generr.data <- function(){
        
        set.seed(2)
        # Define system dimensions
        N <- 50
        sites <- 900
        
        # Environmental predictors for each site
        e1 <- rnorm(sites)
        e2 <- rnorm(sites)
        
        # uncorrelated coefficients for each species and parameters
        sigma1 <- 0.1 # sd beta1
        mean1 <- -1 # mean beta1
        sigma2 <- 0.4 # sd beta2
        mean2 <- 1.5 # mean beta2
        rho <- 0.4 # correlation between betas
        
        # Coefficients for generating the variance-covariance matrix
        nu <- 0.3
        s <- 0.5
        Dis <- (as.matrix(dist(1:N))/N)
        
        # coefficients for each species
        Sigma <- nu*exp(-1/(s*s)*(Dis^2)) + diag(N)*sigma1
        z1 <- mvrnorm(mu = rep(0, times = N), Sigma = Sigma)
        
        # Generate correlations
        beta1 <- z1
        
        # vec <- c(1:round(N*0.5), 1:round(N*0.5))
        # Dis_sigma <- as.matrix(dist(vec))+1-diag(N)
        # Dis_sigma <- (Dis_sigma/max(Dis_sigma))
        Dis_sigma <- (as.matrix(dist(1:N))/N)
        
        Sigma <- 1*exp(-1/(0.3*0.3)*(Dis_sigma^2)) + diag(N)*0.2
        sigma_beta1 <- exp(mvrnorm(mu = rep(0.5, times = N), Sigma = Sigma))
        sigma_test <- sigma_beta1
        # Sigma <- 1*exp(-1/(0.2*0.2)*(Dis^2)) + diag(N)*0.1
        # alpha <- exp(mvrnorm(mu = rep(0, times = N), Sigma = Sigma))
        alpha <- exp(rnorm(N, -0.5,0.5))
        nu <- abs(rnorm(N,0.5, 0.1))+1
        lambda <- inv_logit(rnorm(N, 1.1,0.4))*2-1
        
        sigma_beta1 <- sigma_beta1 * sqrt((pi*(1+3*lambda*lambda)*gamma(3/nu)-(16**(1/nu))*lambda*lambda*gamma(0.5+1/nu)*gamma(0.5+1/nu)*gamma(1/nu))/(pi*gamma(1/nu)))
        beta1 <- beta1 - 2**(2.0/nu)*lambda*gamma(0.5+1.0/nu)/(sqrt(pi)*sigma_beta1)
        
        # Simulate data
        dataset <- expand.grid(site=1:sites, id=1:N)
        dataset$S1 <- e1[dataset$site]
        dataset$S2 <- e2[dataset$site]
        dataset$beta1 <- beta1[dataset$id] 
        dataset$alpha <- alpha[dataset$id]
        dataset$nu <- nu[dataset$id]
        dataset$lambda <- lambda[dataset$id]
        dataset$sigma_beta1 <- sigma_beta1[dataset$id]
        
        dataset$p <- exp(-alpha[dataset$id] - (sigma_beta1[dataset$id]*abs(dataset$S1-beta1[dataset$id])/(1+lambda[dataset$id]*sign(dataset$S1-beta1[dataset$id])))**nu[dataset$id])
        
        dataset$obs <- rbinom(n = length(dataset$S1), size = 1, prob = dataset$p)
        dataset <- data.frame(id=dataset$id, p=dataset$p, obs=dataset$obs, alpha=dataset$alpha, beta1=dataset$beta1, sigma_beta1=dataset$sigma_beta1, nu=dataset$nu, lambda=dataset$lambda, S1=dataset$S1, S2=dataset$S2)                
        return(list(dataset=dataset, corr=Dis, corr2=Dis_sigma, corr3=Dis))
}


# Generate fake data to test the extent to which the model works
simulated.data.categorical <- function(){
        
        # Define system dimensions
        N <- 50
        sites <- 300
        
        # Environmental predictors for each site
        e1 <- rnorm(sites)
        e2 <- rnorm(sites)
        
        # uncorrelated coefficients for each species and parameters
        sigma1 <- 0.3 # sd beta1
        mean1 <- -1 # mean beta1
        sigma2 <- 0.4 # sd beta2
        mean2 <- 1.5 # mean beta2
        rho <- 0.4 # correlation between betas
        
        # Coefficients for generating the variance-covariance matrix
        nu <- 3
        s <- 0.5
        Dis <- (as.matrix(dist(1:N))/N)
        alpha4 <- 1
        alpha3 <- 0.5
        alpha2 <- -0.05
        alpha1 <- -1
        
        # coefficients for each species
        Sigma <- nu*exp(-1/(s*s)*(Dis^2)) + diag(N)*sigma1*sigma1
        z1 <- mvrnorm(mu = rep(mean1, times = N), Sigma = Sigma)
        
        # Generate correlations
        beta1 <- z1
        
        vec <- c(1:round(N*0.5), 1:round(N*0.5))
        Dis_sigma <- as.matrix(dist(vec))+1-diag(N)
        Dis_sigma <- (Dis_sigma/max(Dis_sigma))
        
        Sigma <- 1*exp(-1/(0.3*0.3)*(Dis_sigma^2)) + diag(N)*0.1
        sigma_beta1 <- exp(mvrnorm(mu = rep(0, times = N), Sigma = Sigma))
        Sigma <- 1*exp(-1/(0.2*0.2)*(Dis^2)) + diag(N)*0.1
        # alpha <- exp(mvrnorm(mu = rep(0, times = N), Sigma = Sigma))
        alpha <- exp(rnorm(N, -0.4,1))
        
        # Simulate data
        dataset <- expand.grid(site=1:sites, id=1:N)
        dataset$S1 <- e1[dataset$site]
        dataset$S2 <- e2[dataset$site]
        dataset$beta1 <- beta1[dataset$id] 
        dataset$alpha <- alpha[dataset$id]
        dataset$sigma_beta1 <- sigma_beta1[dataset$id]
        
        dataset$alpha1 <- rep(alpha1, length(dataset$id))
        dataset$alpha2 <- rep(alpha2, length(dataset$id))
        dataset$alpha3 <- rep(alpha3, length(dataset$id))
        dataset$alpha4 <- rep(alpha4, length(dataset$id))
        
        dataset$ppp <- exp(- alpha[dataset$id] - sigma_beta1[dataset$id]*(beta1[dataset$id] - dataset$S1)**2)
        dataset$prob_2to5 <- exp(-exp(alpha1) - alpha[dataset$id] - sigma_beta1[dataset$id]*(beta1[dataset$id] - dataset$S1)**2)              
        dataset$prob_3to5 <- exp(-exp(alpha2) - alpha[dataset$id] - sigma_beta1[dataset$id]*(beta1[dataset$id] - dataset$S1)**2)
        dataset$prob_4to5 <- exp(-exp(alpha3) - alpha[dataset$id] - sigma_beta1[dataset$id]*(beta1[dataset$id] - dataset$S1)**2)
        dataset$prob_5 <- exp(-exp(alpha4) - alpha[dataset$id] - sigma_beta1[dataset$id]*(beta1[dataset$id] - dataset$S1)**2)
        dataset$prob_1 <- 1 - dataset$prob_2to5
        dataset$prob_2 <- dataset$prob_2to5 - dataset$prob_3to5
        dataset$prob_3 <- dataset$prob_3to5 - dataset$prob_4to5
        dataset$prob_4 <- dataset$prob_4to5 - dataset$prob_5
        
        obs <- c()
        for (i in 1:length(dataset$prob_4)) {
                obs[i] <- sample(
                        x = c(1:5), 
                        size = 1, 
                        prob = c(dataset$prob_1[i], dataset$prob_2[i], dataset$prob_3[i], dataset$prob_4[i], dataset$prob_5[i])
                )
        }
        
        dataset$obs <- obs
        dataset <- data.frame(id=dataset$id, obs=dataset$obs, alpha=dataset$alpha, alpha1=dataset$alpha1, alpha2=dataset$alpha2, alpha3=dataset$alpha3, alpha4=dataset$alpha4,  beta1=dataset$beta1, sigma_beta1=dataset$sigma_beta1, S1=dataset$S1, S2=dataset$S2)
        
        # check if all species have at least a couple of occurrences
        test <- colSums(matrix(dataset$obs, sites, sp))
        if(any(test<2)){
                warning('There is at least one species that have less than two occurrences. You might want to rerun (try different parameter values)')
        }
        
        return(list(dataset=dataset, corr=Dis, corr2=Dis_sigma, corr3=Dis))
}

# Generate fake data to test the extent to which the model works
simulated.data.skew <- function(){
        
        # Define system dimensions
        N <- 50
        sites <- 300
        
        # Environmental predictors for each site
        e1 <- rnorm(sites)
        e2 <- rnorm(sites)
        
        # uncorrelated coefficients for each species and parameters
        sigma1 <- 0.3 # sd beta1
        mean1 <- -1 # mean beta1
        sigma2 <- 0.4 # sd beta2
        mean2 <- 1.5 # mean beta2
        rho <- 0.4 # correlation between betas
        
        # Coefficients for generating the variance-covariance matrix
        nu <- 3
        s <- 0.5
        Dis <- (as.matrix(dist(1:N))/N)
        
        # coefficients for each species
        Sigma <- nu*exp(-1/(s*s)*(Dis^2)) + diag(N)*sigma1
        z1 <- mvrnorm(mu = rep(mean1, times = N), Sigma = Sigma)
        
        # Generate correlations
        beta1 <- z1
        
        vec <- c(1:round(N*0.5), 1:round(N*0.5))
        Dis_sigma <- as.matrix(dist(vec))+1-diag(N)
        Dis_sigma <- (Dis_sigma/max(Dis_sigma))
        
        Sigma <- 1*exp(-1/(0.3*0.3)*(Dis_sigma^2)) + diag(N)*0.1
        sigma_beta1 <- exp(mvrnorm(mu = rep(0, times = N), Sigma = Sigma))
        # Sigma <- 1*exp(-1/(0.2*0.2)*(Dis^2)) + diag(N)*0.1
        # alpha <- exp(mvrnorm(mu = rep(0, times = N), Sigma = Sigma))
        alpha <- exp(rnorm(N, -1.5,0.4))
        lambda <- rnorm(N, 0,5)
        # lambda <- rep(5, N)
                
        lambda_hat <- lambda/sqrt(1+lambda**2)
        sigma_hat <- sigma_beta1 * (1 - (2*(lambda_hat**2))/pi)
        beta_hat <- beta1 - sqrt(1/(2*sigma_hat)) * lambda_hat * sqrt(2/pi)
        
        # Simulate data
        dataset <- expand.grid(site=1:sites, id=1:N)
        dataset$S1 <- e1[dataset$site]
        dataset$S2 <- e2[dataset$site]
        dataset$beta1 <- beta1[dataset$id]
        dataset$beta_hat <- beta_hat[dataset$id] 
        dataset$alpha <- alpha[dataset$id]
        dataset$lambda <- lambda[dataset$id]
        
        dataset$sigma_beta1 <- sigma_beta1[dataset$id]
        dataset$sigma_hat <- sigma_hat[dataset$id]
        
        delta <- lambda[dataset$id]/sqrt(1+lambda[dataset$id]**2)
        mu_z <- sqrt(2/pi)*delta
        maxy = 0.5 * ( 4 - pi ) * (delta * sqrt(2/pi))**3 / (1 - 2 * delta**2 / pi )**(3 / 2.0);
        maxy = beta_hat[dataset$id] + 1 / sqrt( 2 * sigma_hat[dataset$id]) * (mu_z - maxy * sqrt(1 - mu_z**2 ) * 0.5 - 0.5 * sign(lambda[dataset$id]) * exp(- 2 * pi / abs(lambda[dataset$id]) ))
        maxy = exp(- sigma_hat[dataset$id] * (maxy - beta_hat[dataset$id])**2) * (1 + pracma::erf((lambda[dataset$id] * (maxy - beta_hat[dataset$id])) * sqrt(sigma_hat[dataset$id]) ))

        dataset$alpha_hat <- log(maxy)+alpha[dataset$id]
        
        dataset$p <- exp(-dataset$alpha_hat - sigma_hat[dataset$id]*(beta_hat[dataset$id] - dataset$S1)**2) * (1 + pracma::erf(lambda[dataset$id] * (dataset$S1-beta_hat[dataset$id]) * sqrt(sigma_hat[dataset$id])))
        
        dataset$obs <- rbinom(n = length(dataset$S1), size = 1, prob = dataset$p)
        dataset <- data.frame(id=dataset$id, obs=dataset$obs, p=dataset$p,
                              alpha=dataset$alpha, alpha_hat=dataset$alpha_hat, beta1=dataset$beta1, beta_hat=dataset$beta_hat,
                              sigma_beta1=dataset$sigma_beta1, sigma_hat=dataset$sigma_hat,
                              lambda = dataset$lambda, S1=dataset$S1, S2=dataset$S2)                
        return(list(dataset=dataset, corr=Dis, corr2=Dis_sigma, corr3=Dis))
}


visua <- function(dataset){
        library(RColorBrewer)
        id = unique(dataset$id)
        n <- length(id)
        qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
        col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
        for(i in 1:n){
                j <- id[i]
                dat <- dataset[dataset$id==j,]
                index <- sort(dat$S1, index.return=T)
                if(i==1){
                        plot(dat$S1[index$ix], dat$p[index$ix], type="l", col=col_vector[j], ylim=c(0,1), ylab="density", xlab="variable")
                }else{
                        lines(dat$S1[index$ix], dat$p[index$ix], type="l", col=col_vector[j])
                }
                abline(v=dat$beta1[1], col=col_vector[j])
        }
        hist(unique(dataset$beta_hat))
        hist(unique(dataset$sigma_beta1))
        hist(unique(dataset$alpha))
}

####
# Run species distribution model for a given species
####
species_distribution.data <- function(variables=c("bio5_", "bio6_","bio12_", "gdd5_", "bio1_","bio15_","bio17_", "bio8_", "TabsY_"),
                                      pca=F, ndim=1,
                                      simulated=F, simulated.type="linear.corr", min.occurrence=0, elevation=F){
        
        if(simulated){
                if(!(simulated.type %in% c("skew", "skew.generror", "generror", "categorical","linear.corr", "linear.gauss", "linear.corr.gauss", "gauss.gauss"))){
                        stop(paste("'", simulated.type, "' is not a valid 'simulated.type'", sep=""))
                }
                if(simulated.type=="categorical"){
                        dataset <- simulated.data.categorical()
                }else if (simulated.type=="skew"){
                        dataset <- simulated.data.skew()
                }else if (simulated.type=="generror"){
                        dataset <- simulated.generr.data()
                }else if (simulated.type=="skew.generror"){
                        dataset <- simulated.skew.generr.data()
                }else{
                        dataset <- simulated.data(simulated.type=simulated.type)
                }
                return(dataset)
        }else{
                if(pca){
                        variables=c("bio5_", "bio6_","bio12_", "gdd5_", "bio1_","bio15_","bio17_", "bio8_", "TabsY_")  
                }
                
                if(elevation){
                        variables=c("bio1_")  
                }
                
                # Load all data
                dat <- prepare.data(variables = variables, min.occurrence=min.occurrence, elevation=elevation)
                
                # extract environmental data
                clim <- as.data.frame(raster::extract(dat$bioclim.data,data.frame(easting = dat$obs.data$easting, northing = dat$obs.data$northing)))
                
                # Find main axes
                if(pca){
                        # colnames(clim) <- unlist(lapply(strsplit(colnames(clim), split = "_"), function(xx) xx[1]))
                        pca.clim <- prcomp(clim, center = TRUE, scale = TRUE) 
                        # Prepare full dataset
                        dataset = cbind(dat$obs.data, pca.clim$x[,c(1:ndim)])
                        
                        # plist <- list()
                        # plist[[1]] <- fviz_eig(pca.clim, barcolor = "gray", barfill = "gray")+
                        #         scale_y_continuous(expand = expansion(add = c(0, 0)))+
                        #         theme_bw()+
                        #         xlab("dimensions")+ylab("percentage of explained variances")+
                        #         theme(plot.title = element_blank(), text = element_text(size=11),
                        #               panel.grid.major = element_blank(),
                        #               panel.grid.minor = element_blank())
                        # p <- fviz_pca_var(pca.clim,
                        #              col.var = "contrib", # Color by contributions to the PC
                        #              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                        #              repel = T,     # Avoid text overlapping
                        # ) +
                        #         xlab("dim 1 (82.6%)")+ylab("dim 2 (9.5%)")+
                        #         theme(plot.title = element_blank(), text = element_text(size=11),
                        #           panel.grid.major = element_blank(),
                        #           panel.grid.minor = element_blank(), legend.title = element_blank())
                        # plist[[2]] <- p+theme(legend.position = "none")
                        # 
                        # grobs <- list()
                        # widths <- list()
                        # heights <- list()
                        # 
                        # for (k in 1:length(plist)){
                        #         grobs[[k]] <- ggplotGrob(plist[[k]])
                        #         widths[[k]] <- grobs[[k]]$widths[2:5]
                        #         heights[[k]] <- grobs[[k]]$heights[2:5]
                        # }
                        # maxwidth <- do.call(grid::unit.pmax, widths)
                        # maxheight <- do.call(grid::unit.pmax, heights)
                        # for (k in 1:length(grobs)){
                        #         grobs[[k]]$widths[2:5] <- as.list(maxwidth)
                        #         grobs[[k]]$heights[2:5] <- as.list(maxheight)
                        # }
                        # 
                        # grobs[[3]] <- get_legend(p)
                        # 
                        # p <- grid.arrange(grobs=grobs, ncol=3, nrow=1, widths=c(0.8,1,0.2))
                        
                }else{
                        if(elevation){
                                # Prepare full dataset
                                dataset = cbind(dat$obs.data, data.frame(elevationstd=dat$obs.data$elevation))
                        }else{
                                # Prepare full dataset
                                dataset = cbind(dat$obs.data, clim)
                        }
                }
                
                # Standarize environmental variables
                for(i in c(1:ncol(dataset))[-c(1:(ncol(dataset)-ndim))]){dataset[,i] <- as.vector(scale(dataset[,i]))}
                
                return(list(dataset=dataset, corr=dat$denv, corr2=dat$dvar, corr3=dat$dtrait))
        }
}

playing.with.multivariate <- function(n, s, c, mu, N){
        # vec <- c(1:round(N*0.5), 1:round(N*0.5))
        vec <- 1:N
        Dis <- dist(vec)
        Sigma <- n*exp(-(1/(s*s))*((as.matrix(Dis)/max(Dis))^2)) + diag(N)*c
        z1 <- mvrnorm(mu = rep(mu, times = N), Sigma = Sigma)
        plot(1:N, z1)
        return(z1)
}


library("sp")
library("raster")
library("maptools")
library("dismo")
library("rJava")
library(SDMTools)
source("../prepare-data/useful-tools.R")
library(MASS)
library(gridExtra)
library(grid)

# Define projection as global variable
projection <- "+proj=somerc +init=world:CH1903"

prepare.data <- function(variables = c("bio5_", "bio6_","bio12_"), min.occurrence=0, elevation=F){
        
        # Determine geographic extent of our data
        places <- read.csv(file = "../../data/properties/codes/places_codes.csv", header = T)
        # places$elevation <- as.vector(scale(places$elevation))
        max.lat <- max(places$northing)
        min.lat <- min(places$northing)
        max.lon <- max(places$easting)
        min.lon <- min(places$easting)
        geographic.extent <- extent(x = c(min.lon, max.lon, min.lat, max.lat))
        
        # Load data for all species
        files <- as.numeric(gsub(".csv", "", list.files("../../data/processed/sdm/")))
        sp_codes <- read.table("../../data/properties/codes/sp_codes.csv", sep=",", header = T)
        dictionary <- read.table("../../data/properties/codes/dictionary.csv", sep=",", header = T)
        correlation_matrix_ids <- read.table("../../data/properties/codes/correlation_matrix_ids.csv", sep="\t", header = F)
        denvironment <- as.matrix(read.table("../../data/properties/distance-matrices/environment.csv", sep=","))
        dvariation <- as.matrix(read.table("../../data/properties/distance-matrices/variation.csv", sep=","))
        dtraits <- as.matrix(read.table("../../data/properties/distance-matrices/trait.csv", sep=","))
        
        indicator <- read.table("../../data/properties/codes/temperature_indicator.csv", sep=",")
        neophytes <- read.table("../../data/properties/codes/neophytes-list.csv", sep=",")
        tendency <- read.table("../../data/properties/codes/change-tendency.csv", sep=",")
        competitive <- read.table("../../data/properties/codes/competitive_indicator.csv", sep=",")
        detailed <- read.table("../../data/properties/codes/neophytes-detailed.csv", sep=",")
        
        # Check that there aren't unnexpected files
        if(!all(sort(files)==1:length(files))){
            stop("Odd files in the folder")    
        }
        
        # Prepare main file
        obs.data <- data.frame()
        name.idx <- c()
        kdx <- 1
        Tind <- c()
        NEO <- c()
        Tend <- c()
        compet <- c()
        deta <- c()
        
        # Read observations
        for(idx in 1:length(files)){
             if(sp_codes$range[sp_codes$id==idx]<min.occurrence){
                     next
             }
             obs.data_ <- read.csv(file = paste("../../data/processed/sdm/", as.character(idx), ".csv", sep = ""), header = T)
             ### Ok so I think the order for the correlation matrix is the same, but I do need to double-check
             ### This is a very dumb way of doing just that. I didn't want to think.
             new.name <- as.character(dictionary$new.names[as.character(dictionary$old.names)==as.character(sp_codes$sp[idx])])
             Tind <- rbind(Tind, c(kdx, new.name, as.character(indicator$nflor.T[new.name==indicator$nflor.spnames])))
             NEO <- rbind(NEO, c(kdx, new.name, neophytes$neo[new.name==neophytes$names]))
             Tend <- rbind(Tend, c(kdx, new.name, tendency$decrease[new.name==tendency$names], tendency$decrease.low[new.name==tendency$names], tendency$increase[new.name==tendency$names], tendency$other[new.name==tendency$names]))
             compet <- rbind(compet, c(kdx, new.name, as.character(competitive$nflor.KS[new.name==competitive$nflor.spnames])))
             deta <- rbind(deta, c(kdx, new.name, detailed$i[new.name==detailed$names],
                                   detailed$n[new.name==detailed$names],
                                   detailed$nw[new.name==detailed$names],
                                   detailed$a[new.name==detailed$names],
                                   detailed$ja[new.name==detailed$names], 
                                   detailed$jn[new.name==detailed$names],
                                   detailed$isn[new.name==detailed$names],
                                   detailed$isa[new.name==detailed$names],
                                   detailed$isj[new.name==detailed$names],
                                   detailed$isasn[new.name==detailed$names],
                                   detailed$asn[new.name==detailed$names]))
             
             name.idx <- c(name.idx,correlation_matrix_ids$V1[as.character(correlation_matrix_ids$V2)==new.name])
             ###
             obs.data_$id <- kdx
             obs.data_$real.id <- idx
             obs.data_$mmsbm.id <- correlation_matrix_ids$V1[as.character(correlation_matrix_ids$V2)==new.name]
             if(elevation){
                     obs.data_$elevation <- places$elevation
             }
             obs.data <- rbind(obs.data, obs.data_)
             kdx <- kdx+1
        }
        
        newdeta <- as.data.frame(deta)
        colnames(newdeta) <- c("id", colnames(detailed))
        # class(deta) <- "numeric"
        # newdeta <- newdeta[,c(T, T, colSums(deta[,3:ncol(deta)])>0)]

        
        if(min.occurrence==10){
                write.table(Tind, "../../data/properties/codes/temperature_indicator_reindexed.csv", sep=",")
                write.table(NEO, "../../data/properties/codes/neophytes-list_reindexed.csv", sep=",")
                write.table(Tend, "../../data/properties/codes/change-tendency_reindexed.csv", sep=",")
                write.table(compet, "../../data/properties/codes/competitive_reindexed.csv", sep=",")
                write.table(newdeta, "../../data/properties/codes/neophytes-detailed_reindexed.csv", sep=",")
                
        }else{
                write.table(Tind, paste("../../data/properties/codes/temperature_indicator_reindexed-",as.character(min.occurrence),".csv", sep=""))
                write.table(NEO, paste("../../data/properties/codes/neophytes-list_reindexed-",as.character(min.occurrence),".csv", sep=""))
                write.table(Tend, paste("../../data/properties/codes/change-tendency_reindexed-",as.character(min.occurrence),".csv", sep=""))
                write.table(compet, paste("../../data/properties/codes/competitive_reindexed-",as.character(min.occurrence),".csv", sep=""))
                write.table(newdeta, paste("../../data/properties/codes/neophytes-detailed_reindexed-",as.character(min.occurrence),".csv", sep=""))
                
        }
        
        # reshape correlation matrices
        denvironment <- denvironment[,name.idx]
        denvironment <- denvironment[name.idx,]
        dvariation <- dvariation[,name.idx]
        dvariation <- dvariation[name.idx,]
        dtraits <- dtraits[,name.idx]
        dtraits <- dtraits[name.idx,]
        
        # rename cols and rows
        colnames(denvironment) <- 1:ncol(denvironment)
        rownames(denvironment) <- 1:nrow(denvironment)
        colnames(dvariation) <- 1:ncol(dvariation)
        rownames(dvariation) <- 1:nrow(dvariation)
        colnames(dtraits) <- 1:ncol(dtraits)
        rownames(dtraits) <- 1:nrow(dtraits)
        
        obs.data$obs <- 1*(obs.data$abundance>0)
        obs.data$abundance <- factor(obs.data$abundance,
                               levels = sort(unique(obs.data$abundance)),
                               labels = 1:length(unique(obs.data$abundance)))

        # Load environmental data
        files <- list.files(path="../../data/raw/climatic-data/", pattern = "bil$", full.names = TRUE)
        files <- files[grepl(paste(variables,collapse="|"), files)]
        bioclim.data <- stack(files)
        crs(bioclim.data)<-projection
        
        return(list(obs.data = obs.data, bioclim.data = bioclim.data, xlim = c(min.lon, max.lon), ylim = c(min.lat, max.lat), 
                    denv = denvironment, dvar = dvariation, dtrait=dtraits))
}

# Generate fake data to test the extent to which the model works
simulated.data <- function(simulated.type="linear.corr"){
        
        set.seed(3)
        
        # Define system dimensions
        N <- 25
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

        # vec <- c(1:round(N*0.5), 1:round(N*0.5))
        # Dis_sigma <- as.matrix(dist(vec))+1-diag(N)
        # Dis_sigma <- (Dis_sigma/max(Dis_sigma))
        Dis_sigma <- (as.matrix(dist(1:N))/N)
        
        Sigma <- 1*exp(-1/(0.3*0.3)*(Dis_sigma^2)) + diag(N)*0.1
        sigma_beta1 <- exp(mvrnorm(mu = rep(0, times = N), Sigma = Sigma))
        # Sigma <- 1*exp(-1/(0.2*0.2)*(Dis^2)) + diag(N)*0.1
        # alpha <- exp(mvrnorm(mu = rep(0, times = N), Sigma = Sigma))
        alpha <- exp(rnorm(N, -0.5,1))

        # Simulate data
        dataset <- expand.grid(site=1:sites, id=1:N)
        dataset$S1 <- e1[dataset$site]
        dataset$S2 <- e2[dataset$site]
        dataset$beta1 <- beta1[dataset$id] 
        dataset$alpha <- alpha[dataset$id]
        dataset$sigma_beta1 <- sigma_beta1[dataset$id]
        
        dataset$p <- exp(-alpha[dataset$id] - sigma_beta1[dataset$id]*(beta1[dataset$id] - dataset$S1)**2)
        
        dataset$obs <- rbinom(n = length(dataset$S1), size = 1, prob = dataset$p)
        dataset <- data.frame(id=dataset$id, obs=dataset$obs, alpha=dataset$alpha, beta1=dataset$beta1, sigma_beta1=dataset$sigma_beta1, S1=dataset$S1, S2=dataset$S2)                
        return(list(dataset=dataset, corr=Dis, corr2=Dis_sigma, corr3=Dis))
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


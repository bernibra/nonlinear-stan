# Simulate data for a skewed generalized error distribution. Function used to create the first figure of the manuscript. 
generrskew <- function(x, alpha, beta, gamma, lambda, nu){

  # Transform gamma and beta
  gamma_v <- sqrt(gamma * ((pi*(1+3*lambda*lambda)*base::gamma(3/nu)-(16**(1/nu))*lambda*lambda*base::gamma(0.5+1/nu)*base::gamma(0.5+1/nu)*base::gamma(1/nu))/(pi*base::gamma(1/nu))))
  beta_m <- beta - 2**(2.0/nu)*lambda*base::gamma(0.5+1.0/nu)/(sqrt(pi)*gamma_v)
  
  return(exp(-alpha - (gamma_v*abs(x-beta_m)/(1+lambda*sign(x-beta_m)))**nu))
}

# Calculate the excess kurtosis for a skewed generalized error distribution.
kurtosis.skew.generror <- function(nu, lambda){
  v <- ((pi*(1+3*lambda*lambda)*gamma(3/nu)-((16)**(1/nu))*lambda*lambda*(gamma(0.5+1/nu)**2)*gamma(1/nu))/(pi*gamma(1/nu)))**(-1/2)
  kurtosis <- ((v**4)/(pi*pi*gamma(1/nu)))*(-3*(256**(1/nu))*(lambda**4)*(gamma(0.5+1/nu)**4)*gamma(1/nu) +
                                              3*(2**((4+nu)/nu))*pi*(lambda**2)*(1+3*(lambda**2))*(gamma(0.5+1/nu)**2)*gamma(3/nu)-
                                              (2**(4+2/nu))*(pi**(3/2))*(lambda**2)*(1+lambda**2)*gamma(0.5+1/nu)*gamma(4/nu)+
                                              (pi**2)*(1+10*(lambda**2)+5*(lambda**4))*gamma(5/nu))-3
  return(kurtosis)
}

# Calculate the skewness for a skewed generalized error distribution.
skewness.skew.generror <- function(nu, lambda){
  v <- ((pi*(1+3*lambda*lambda)*gamma(3/nu)-((16)**(1/nu))*lambda*lambda*(gamma(0.5+1/nu)**2)*gamma(1/nu))/(pi*gamma(1/nu)))**(-1/2)
  skewness <- (lambda*(v**3)/((pi**(3/2))*gamma(1/nu)))*((2**((6+nu)/nu))*(lambda**2)*(gamma(0.5+1/nu)**3)*gamma(1/nu)-
                                                           3*(4**(1/nu))*pi*(1+3*lambda*lambda)*gamma(0.5+1/nu)*gamma(3/nu)+
                                                           4*(pi**(3/2))*(1+lambda*lambda)*gamma(4/nu))
  return(skewness)
}

# Code to generate figure 1
plot.one.distribution <- function(alpha, beta, gamma, lambda, nu){
  x <- seq(-4, 4, length.out = 2000)
  dat <- data.frame(x=x, y=generrskew(x, alpha, beta, gamma, lambda, nu))
  
  return(ggplot(data = dat, aes(x=x, y=y)) +
           geom_line(color="black")+
           theme_bw() +
           xlab("x")+
           ylab("probability")+
           coord_cartesian(ylim=c(0, 1), xlim=c(-4, 4)) +
           theme(text = element_text(size=10),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 legend.position="none",
                 panel.border = element_rect(colour = "black", fill=NA, size=0.3),
                 plot.title = element_text(size=10)))
}

# Code to generate figure 2
plot.ranking.x <- function(mu, ci, color="#d95f02"){
  # Sort data
  mu_order <- sort(mu,index.return=T)$ix

  # Generate data.frame  
  df <- data.frame(sp=1:length(mu), mean=mu[mu_order], low=ci[1,mu_order], high=ci[2,mu_order])
  
  p <- ggplot(df, aes(x=sp, y=mean)) + 
    geom_segment(aes(x=sp, xend=sp, y=low, yend=high), data=df, color=color, alpha=0.5) +
    geom_point(color=color, alpha=0.7)+
    ylab("values")+
    xlab("species")+
    scale_x_continuous(expand = expansion(add = c(0.2, 0.2)), breaks= pretty_breaks())+
    scale_y_continuous(expand = expansion(add = c(0, 0)))

  return(p)
}

log.vs.probability <- function(X, mfit, samples = 2000){
  loglik <- mfit$draws(c("log_lik"), format = "data.frame") %>% select(-c(`.draw`, `.chain`, `.iteration`))
  gammav <- mfit$draws(c("gammav"), format = "data.frame") %>% select(-c(`.draw`, `.chain`, `.iteration`))
  betam <- mfit$draws(c("betam"), format = "data.frame") %>% select(-c(`.draw`, `.chain`, `.iteration`))
  lambda <- mfit$draws(c("lambda"), format = "data.frame") %>% select(-c(`.draw`, `.chain`, `.iteration`))
  nu <- mfit$draws(c("nu"), format = "data.frame") %>% select(-c(`.draw`, `.chain`, `.iteration`))
  
  maxlog <- max(loglik)
  minlog <- min(loglik)
  
  dat <- expand.grid(seq(from=0, to=1, length.out = 100),
                     seq(from=minlog, to=maxlog, length.out = 100))
  
  z.hex <- hexbin(dat$Var1, dat$Var2, xbnds = c(0, 1), ybnds = c(minlog, maxlog), xbins = 25, IDs = T)
  data <- rep(0,max(as.numeric(names(table(z.hex@cID)))))
  finaldata <- data.frame(hexbin::hcell2xy(z.hex),
                          cell = z.hex@cell,
                          count = z.hex@count)
  normalization <- finaldata
  
  for(i in 1:samples){
    loglik_ <- loglik[i,]
    betam_ <- outer(X, as.matrix(betam[i,]), "-")
    gammav_ <- outer(rep(1, length(X)),as.matrix(gammav[i,]), "*")
    lambda_ <- outer(rep(1, length(X)),as.matrix(lambda[i,]), "*")
    nu_ <- outer(rep(1, length(X)),as.matrix(nu[i,]), "*")
    
    p <- as.vector(exp(-(gammav_ * abs(betam_)/(1+lambda_*sign(betam_)))**(nu_)))
    
    z.hex_ <- hexbin(p, loglik_, xbnds = c(0, 1), ybnds = c(minlog, maxlog), xbins = 25, IDs = T)
    z.hex_ <- table(z.hex_@cID)
    z.hex_id <- as.numeric(names(z.hex_))
    
    for(k in 1:length(z.hex_)){
      data[z.hex_id[k]] <- data[z.hex_id[k]] + z.hex_[k]
    }
    
  }
  
  finaldata$count <- data
  finaldata <- finaldata[finaldata$count!=0,]
  finaldata <- finaldata[finaldata$y<=log(0.50),]
  
  finaldata$count <- finaldata$count/sum(  finaldata$count )
  
  return(finaldata)
}


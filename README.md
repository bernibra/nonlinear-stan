# nonlinear-stan [![DOI](https://zenodo.org/badge/388764760.svg)](https://zenodo.org/badge/latestdoi/388764760) [![](https://badgen.net/badge/DOI/10.1086%2F727518/green)](https://doi.org/10.1086/727518)

Bayesian hierarchical non-linear models to study the shape of plant distributions. If you want to run the models with the data used in the manuscript, this is publicly available [here](https://doi.org/10.1111/j.1472-4642.2011.00792.x).

## How to cite?

Uncovering Broad Macroecological Patterns by Comparing the Shape of Species’ Distributions along Environmental Gradients<br>
Bernat Bramon Mora, Antoine Guisan, and Jake M. Alexander<br>
The American Naturalist 2024 203:1, 124-138

## Considerations

I would not recommend running these models unless you have access to a HPC cluster, as these can be computationally very expensive. All model runs were done using the computer cluster at ETHZ (Euler). This is an estimate of the hardware requirements to run the example below:

Resource summary:
- CPU time :                                   263403.94 sec.
- Max Memory :                                 11969 MB
- Average Memory :                             1251.61 MB
- Max Processes :                              6
- Max Threads :                                49
- Run time :                                   8207 sec.
- Turnaround time :                            166988 sec.

## Download the code

Download or clone this repository in your computer using `git clone https://github.com/bernibra/nonlinear-stan.git`

## Process data

To simulate data, one can use the functions provided in `./code/simulated-data/simulate-data.R`. For example, to generate presence/absence data for 50 species accross 900 sites presenting fat-tailed and skewed distributions, you can run the following:
```r
source("./code/simulate-data.R")

L <- 50
N <- 900

my_data <- simulated.generror(sp=L, sites=N)
```
Notice that this will simulate data according to the parameter default values. You can change any of those following the parameter naming (see the manuscript methods section):
  - `beta_s`: $\sigma_{\beta}$
  - `beta_mu`: $\hat{\beta}$
  - `gamma_s`: $\sigma_{\gamma}$
  - `gamma_mu`: $\hat{\gamma}$
  - `beta_eta`: $\eta_{\beta}$
  - `beta_rho`: $\rho_{\beta}$
  - `gamma_nu`: $\eta_{\gamma}$
  - `gamma_rho`: $\rho_{\gamma}$
  - `alpha_mu`: $\hat{\alpha}$
  - `alpha_s`: $\sigma_{\alpha}$
  - `nu_mu`: $\hat{\nu}$
  - `nu_s`: $\sigma_{\nu}$
  - `lambda_mu`: $\hat{\lambda}$
  - `lambda_s`: $\sigma_{\lambda}$

You can have a quick look at what the theoretical distributions looks like using:
```r
plot.distributions(my_data$dataset)
```

If you are instead using the empirical data used in the manuscript, you can simply create a list with the different elements in the object `my_data`. Fundamentally, this should include a `dataset` data.frame with $sites \times species$ rows, with columns `id` (species integer ID), `obs` (presence/absence data or ordered categorical variable), `S1` (the environmental variable; scaled for optimal sampling), and `site` (site integer ID). Finally, the `my_data` list should also include two distance matrices `Dgamma` and `Dbeta` for the Gaussian process (see manuscript methods section for details).

## Prepare the data for the model

Before sampling the posterior distributions, we need to prepare the data and define the initial values. Notice that we generate two objects, one with the input data and another with the initial values for each chain (here we use 3 chains).

```r
# Load the cmdstanr library
library(cmdstanr)

# Reformat the data
obs <- t(matrix(my_data$dataset$obs, N, L))
X <- matrix(my_data$dataset$S1, N, L)[,1]
indices <- 1:L

# Input data object for stan
dat_1.0 <- list(N=N,
                L=L,
                minp=1e-100,
                Y=obs,
                indices=indices,
                X1=X,
                Dmat_b=my_data$Dbeta,
                Dmat_g=my_data$Dgamma
)

# Set starting values for the parameters
start_1.0 <- list(
  zalpha = rep(0, dat_1.0$L),
  zbeta = rep(0, dat_1.0$L),
  zgamma = rep(0, dat_1.0$L),
  zlambda = rep(0, dat_1.0$L),
  znu = rep(0, dat_1.0$L),
  alpha_bar = 0,
  nu_bar = 0,
  lambda_bar = 0,
  beta_bar = 0,
  gamma_bar = 0,
  sigma_a = 0.1,
  sigma_n = 0.1,
  sigma_l = 0.1,
  sigma_b = 0.1,
  etasq_b = 0.1,
  rhosq_b = 0.1,
  sigma_g = 0.1,
  etasq_g = 0.1,
  rhosq_g = 0.1
)

# Initialize data structure
n_chains_1.0 <- 3
init_1.0 <- list()
for ( i in 1:n_chains_1.0 ) init_1.0[[i]] <- start_1.0
```

## Sample the posterior distributions

The models used in the manuscript can be found in `./code/stan-code/models-binomial.R` and `./code/stan-code/models-categorical.R`. To sample the posterior distributions, we used the R package _cmdstanr_ (see Session Information section at the end of the page). This is because the models use the [`reduce_sum`](https://mc-stan.org/users/documentation/case-studies/reduce_sum_tutorial.html) functionality, a way to parallelise the execution of a single Stan chain across multiple cores. Let's try to do so with the fat-tailed and skewed model.

```r
# Load the stan code
source('./code/stan-code/models-binomial.R')
model_code = binomial.fat.tailed.skewed

# Prepare the model
model <- cmdstan_model(write_stan_file(model_code), cpp_options = list(stan_threads = TRUE))

# Sample the posterior distributions
mfit_1.0 <- model$sample(data = dat_1.0,
                              init = init_1.0,
                              chains = 3,
                              threads_per_chain = 15,
                              parallel_chains = 3,
                              max_treedepth = 15,
                              max_depth = 15,
                              iter_sampling = 1000,
                              refresh = 500)

# Save the object as an rds file
mfit_1.0$save_object(file = "model-cdmstan.rds")
```

Notice that the computational resources requested here are not those of a regular laptop; again, this should be run in a computer cluster where one have a substantial amount of memory and cores. Run it at your own risk on your laptop, lowering the `threads_per_chain` to fit your laptop specifications (usually 1 or 2 per core). If so, it would also be a good idea to lower the number of species `L` and sites `N`.

## Visualize the results

Now, the results of the sampling process are stored in the object `mfit_1.0`, also stored as an RDS object created in the main directory.

### Run extract samples

To visualize the results and produce the figures presented in the manuscript, one needs to first extract the samples. To do so, we will use functions from the R packages `dplyr`, `gridExtra`, `scales`, `rethinking`, `hexbin` and `posterior` (see Session Information section at the end of the page):

```r
library(dplyr)
library(posterior)
library(scales)
library(rethinking)
library(hexbin)
library(gridExtra)

beta <- mfit_1.0$draws(c("beta"), format = "data.frame") %>% select(-c(`.draw`, `.chain`, `.iteration`))
gamma <- mfit_1.0$draws(c("gamma"), format = "data.frame") %>% select(-c(`.draw`, `.chain`, `.iteration`))
alpha <- mfit_1.0$draws(c("alpha"), format = "data.frame") %>% select(-c(`.draw`, `.chain`, `.iteration`))
nu <- mfit_1.0$draws(c("nu"), format = "data.frame") %>% select(-c(`.draw`, `.chain`, `.iteration`))
lambda <- mfit_1.0$draws(c("lambda"), format = "data.frame") %>% select(-c(`.draw`, `.chain`, `.iteration`))
```

### Fig 1
The first figure of the manuscript was generated without data, and it can be reproduced using the function `plot.one.distribution` in `utility.R`

```r
library(ggplot2)
source('./code/utility.R')

plot.one.distribution(alpha = 0.5, beta = 1, gamma = 2, nu = 2.2, lambda = -0.5)
```

### Fig 2
The second figure of the manuscript displays the center and confidence interval of the mean and variance of distributions for each species. We can reproduce this plot for the mean and variance using `plot.ranking.x` in `utility.R`. For example, for the variance:

```r
variance <- 1/gamma

mu_variance <- apply(variance,2,mean)
ci_variance <- apply(variance,2,PI)

g <- plot.ranking.x(mu_variance, ci_variance) +
    coord_trans(y="log", clip = "off")
```

if you are using the simulated data and want to visualize the true values, you can:

```r
mu_order <- sort(mu_variance,index.return=T)$ix

true_variance <- my_data$dataset %>%
                 select(gamma, id) %>%
                 distinct() %>%
                 mutate(variance = 1.0/gamma) %>%
                 pull(variance)

dat <- data.frame(x = 1:length(true_variance), y = true_variance[mu_order])

g + geom_point(data = dat, aes(x=x, y=y), color = "black", shape = 1)
```

### Fig 3
The third figure of the manuscript displays the posterior distribution of the average kurtosis and average skewness across species. One can calculate those using `kurtosis.skew.generror` and `skewness.skew.generror` from `utility.R`:

```r
skewness <- skewness.skew.generror(nu, lambda)
skewness_bar <- apply(skewness, 1, mean)

kurtosis <- kurtosis.skew.generror(nu, lambda)
kurtosis_bar <- apply(kurtosis, 1, mean)

dat <- data.frame(y = skewness_bar, x = kurtosis_bar)

p1 <- ggplot(dat, aes(x=x)) +
    geom_vline(xintercept = 0, color = "#999999", linetype="dashed") +
    stat_density(geom = "line", position = "identity", colour = "#E69F00") +
    xlab("kurtosis") + 
    ggtitle("(a)  average kurtosis of distributions") +
    theme_bw()
    
p2 <- ggplot(dat, aes(x=y)) +
    geom_vline(xintercept = 0, color = "#999999", linetype="dashed") +
    stat_density(geom = "line", position = "identity", color = "#E69F00") +
    xlab("skewness") + 
    ggtitle("(b)  average skewness of distributions") +
    theme_bw()
    
grid.arrange(p1, p2, ncol=2)
```

### Fig 4
The fourth figure studies the correlation across pairs of parameters. To do so, we just need to calculate the mean of the posterior distribution for two given parameters across species, and display them against each other.

```r
mu_beta <- apply(beta,2,mean)

dat <- data.frame(mean=mu_beta, variance=mu_variance)
ggplot(dat, aes(y=mean, x=variance)) + 
  geom_point(color="#1b9e77")+
  coord_trans(x="log", clip = "off") +
  scale_x_continuous(expand = expansion(add = c(0, 0))) +
  scale_y_continuous(expand = expansion(add = c(0, 0))) +
  theme_bw()
```

Notice that for the upper panels and side margin panels comparing groups of species (only one in this case), one can obtain the distribution of means across groups by simply applying:

```r
average_beta <- apply(beta,1,mean)
average_variance <- apply(variance,1,mean)
```

### Fig 5

The final figure studies the log-likelihood and how the errors are distributed across distributions. To reproduce this figure we need to extract the log-likelihood values for every sample and calculate the normalized probability for the distribution of every species. We can do this with the function `log.vs.probability` from `utility.R`:

```r
data <- log.vs.probability(X, mfit_1.0, samples = 2000)

ggplot(data, aes(x=x, y=y, fill=count)) + 
    geom_hex(stat="identity",size=0.2, alpha=1) +
    xlab("normalized probability") +
    ylab("log-likelihood") +
    annotate("text", x = c(0.025, 0.985), y = c(-0.17, -0.17), label = c("tails", "peak"), colour = "#525252", size=3) +
    coord_cartesian(clip = 'off') +
    scale_x_continuous(expand = c(0.02, 0.02)) +
    scale_y_continuous(expand = c(0.02, 0.02), limits = c(min(data$y)-0.2, 0)) +
    scale_fill_distiller(limits= c(-0.002,max(data$count)), palette = "PuRd", direction=1, na.value = "white", labels = as.vector(sapply(c(1, 1+1:5*2), function(ll) paste(ll, "%", sep=""))) , breaks = c(1, 1+1:5*2)*0.01) +
    theme_bw() +
    guides(color = FALSE, size = FALSE) +
    theme(text = element_text(size=11),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.background = element_blank(),
          legend.title = element_blank())
```

## Session Information

```R
R version 4.0.4 (2021-02-15)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Big Sur 10.16

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] gridExtra_2.3        hexbin_1.28.2        scales_1.2.0        
 [4] posterior_1.2.1      dplyr_1.0.8          cmdstanr_0.3.0.9000 
 [7] rethinking_2.13      rstan_2.21.3         StanHeaders_2.21.0-7
[10] ggplot2_3.3.6        MASS_7.3-56         

loaded via a namespace (and not attached):
 [1] shape_1.4.6          tidyselect_1.1.2     xfun_0.30           
 [4] purrr_0.3.4          lattice_0.20-45      colorspace_2.0-3    
 [7] vctrs_0.4.0          generics_0.1.2       stats4_4.0.4        
[10] loo_2.5.1            utf8_1.2.2           rlang_1.0.2         
[13] pkgbuild_1.3.1       pillar_1.7.0         glue_1.6.2          
[16] withr_2.5.0          DBI_1.1.2            distributional_0.3.0
[19] matrixStats_0.61.0   lifecycle_1.0.1      munsell_0.5.0       
[22] gtable_0.3.0         mvtnorm_1.1-3        codetools_0.2-18    
[25] coda_0.19-4          inline_0.3.19        knitr_1.39          
[28] callr_3.7.0          ps_1.6.0             fansi_1.0.3         
[31] Rcpp_1.0.8.3         backports_1.4.1      checkmate_2.0.0     
[34] RcppParallel_5.1.5   jsonlite_1.8.0       abind_1.4-5         
[37] farver_2.1.0         tensorA_0.36.2       processx_3.5.3      
[40] grid_4.0.4           cli_3.2.0            tools_4.0.4         
[43] magrittr_2.0.3       tibble_3.1.6         crayon_1.5.1        
[46] pkgconfig_2.0.3      ellipsis_0.3.2       prettyunits_1.1.1   
[49] assertthat_0.2.1     R6_2.5.1             compiler_4.0.4 
```

The only package that is not currently on CRAN is Richard McElreath’s [`rethinking`](https://github.com/rmcelreath/rethinking) package. This can be installed with `devtools` as follows:

```R
devtools::install_github("rmcelreath/rethinking")
```

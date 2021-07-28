# nonlinear-stan

Bayesian hierarchical non-linear models to study the shape of plant distributions

## Considerations

I would not recommend running these models unless you have access to a HPC cluster, as these can be computationally very expensive. All model runs were done using the computer cluster at ETHZ (Euler). An estimate of the hardware requirements is provided below.

## Download the code

Download or clone this repository in your computer using `git clone https://github.com/bernibra/nonlinear-stan.git`

## Simulate data

To simulate data, one can use the functions provided in `./code/simulated-data/simulate-data.R`. For example, to generate presence/absence data for 50 species accross 900 sites presenting fat-tailed and skewed distributions, you can run the following:
```r
source("./code/simulated-data/simulate-data.R")

L <- 50
N <- 900

my_data <- simulated.generror(sp=L, sites=N)
```

## Prepare the data

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
                Dmat_b=d$Dbeta,
                Dmat_g=d$Dgamma
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

The models used in the manuscript can be found in `./code/stan-code/models-binomial.R` and `./code/stan-code/models-categorical.R`. To sample the posterior distributions, we used the R package _cmdstanr_. This is because the models use the [`reduce_sum`](https://mc-stan.org/users/documentation/case-studies/reduce_sum_tutorial.html) functionality, a way to parallelize the execution of a single Stan chain across multiple cores.

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
mfit_1.0$save_object(file = "model-cdmstan.rds", sep="")
```
Resource summary:
- CPU time :                                   263403.94 sec.
- Max Memory :                                 11969 MB
- Average Memory :                             1251.61 MB
- Max Processes :                              6
- Max Threads :                                49
- Run time :                                   8207 sec.
- Turnaround time :                            166988 sec.

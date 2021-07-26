# nonlinear-stan

Bayesian hierarchical non-linear models to study the shape of plant distributions

## Considerations

I would not recommend running these models unless you have access to a HPC cluster, as these can be computationally very expensive. All model runs were done using the computer cluster from ETHZ (Euler). An estimate of the hardware requirements is provided below.

## Download the code

Download or clone this repository in your computer using git clone https://github.com/bernibra/nonlinear-stan.git

## Simulate data

To simulate data, one can use the functions provided in './code/simulated-data/simulate-data.R'. For example, to generate presence/absence data for 50 species accross 900 sites presenting fat-tailed and skewed distributions, you can run the following:
```r
source("./code/simulated-data/simulate-data.R")

d <- simulated.generror(sp=50, sites=900)
```




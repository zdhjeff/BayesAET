# Bayesian-Adaptive-Enrichment-Design
An R package for Bayesian adaptive enrichment trials with multiple treatment arms and multiple (categorical) subpopulations from the paper "Bayesian adaptive enrichment design for multi-arm clinical trials – an R package"  <br>
Run the following command to install the package and other packages:

```r
library(devtools)
install_github("zdhjeff/BayesAET")
#install the following packages if needed
library(abind)
library(reshape2)
library(rjags)
library(doParallel)
library(MASS)
library(fastDummies)
```
# Example:

```r
library(abind)
library(reshape2)
library(rjags)
library(doParallel)
library(MASS)
library(BayesAET)
library(fastDummies)

#' @param nt number of treatment arms
#' @param ns number of subgroups
#' @param maxN the maximum sample size, trial will stop when achieving this number
#' @param ss.interim.es A list containing the g vectors, with gth vector indicating the accumulated number of sample sizes for each interim analysis for gth subpopulation
#' @param response.type either 'binary'(probability), 'count'(lambda) or 'gaussian'
#' @param mean.response matrix of total effect sizes: a nt (row) * ns (col) matrix, with row number indicating treatment and col number indicating subgroup
#' @param prob.subpopulation the probability of a subject coming from one specific subpopulation. default: rep (1/ns, ns)
#' @param prob.trtarm  the (initial) probability of a subject being assigned to a specific treatment arm. default: rep (1/nt, nt)
#' @param upper upper probability threshold to claim superiority (in a specific subgroup); the treatment arm wins if it goes above this threshold
#' @param lower lower probability threshold to claim superiority (in a specific subgroup); the treatment arm will be dropped if it goes below this threshold
#' @param RAR whether using responsive adaptive randomization (RAR)
#' @param MOR the minimum meaningful outcome threshold for each subgroup
#' @param prob.MOR the probability threshold of being larger than the MOR, treatment arms below this threshold will be dropped
#' @param N Number of MCMC samples
#' @param prior.cov the prior covariance for a multivariate normal prior
#' @param prior.mean the prior mean for a multivariate normal prior

## Gaussian outcome
nt=3
ns=2
BAET.sim(nt, ns,
        ss.interim.es = list(c(30, 50, 100), c(20, 60, 100)),
        response.type = "gaussian",
        sig.e = 10,
        mean.response = matrix(c(seq(nt*ns) ), nrow = nt, ncol = ns, byrow = F),
        prob.subpopulation = rep (1/ns, ns),
        prob.trtarm = rep (1/nt, nt),
        maxN = 100,
        upper = rep (1, ns),
        lower = rep (-0.1, ns),
        rar = T,
        N.MCMC = 3000,
        prior.cov = diag(100, ns*nt),
        prior.mean = rep(0, ns*nt))
 

#' @return n.interim: number of interim looks conducted to end the whole trial
#' @return trt_sub: simulated treatment arm allocation (1st column) and subgroup (2nd column)
#' @return est: treatment effect estimation(posterior mean) and the 95% Credible interval bounds for each treatment arm in each subgroup
#' @return powerind: power indicator of whether the best treatment arms is correctly selected in each subgroup
#' @return y: simlulated outcome
#' @return N_terminate: the total sample size consumed when trial ends
#' @return prob_sup_minioutcome: the probability of a treatment large than the minioutcome
#' @return prob_superiority: the probability of a treatment being the best among each subgroup

## outputs:

$n.analysis
[1] 4

$trt_sub # only shows the first 5 samples
trtarm_ind_b subpop_ind_b
[1,]            3            1
[2,]            2            1
[3,]            2            2
[4,]            2            2
[5,]            1            1



$est
$est[[1]]
trt.est  lowbound  upbound   sd.est
coef_all[1]  0.3900897 -3.880923 5.032788 2.272657
coef_all[2] -3.3367884 -8.667656 1.680998 2.646658
coef_all[3]  5.7576970  2.158152 9.305612 1.842488

$est[[2]]
trt.est  lowbound   upbound   sd.est
coef_all[4] 3.036873 -2.679267  8.574124 2.918975
coef_all[5] 6.009781  1.837800  9.895292 2.040647
coef_all[6] 7.026762  2.830796 11.662915 2.244942


$powerind
[1] 0 0

$y # only shows the first 5 samples
[1] 5.29427406 2.89057932 1.05350922 27.27202256 -3.81152040

$N_terminate
[1] 100

$ss.sub
[,1]
[1,]   54
[2,]   46

$prob_sup_minioutcome
$prob_sup_minioutcome[[1]]
[,1]      [,2]      [,3]
[1,] 0.5316667 0.5503333 0.5600000
[2,] 0.2876667 0.1416667 0.1073333
[3,] 0.9460000 0.9910000 1.0000000

$prob_sup_minioutcome[[2]]
[,1]      [,2]
[1,] 0.4433333 0.8496667
[2,] 0.8020000 0.9976667
[3,] 0.8333333 1.0000000


$prob_superiority
$prob_superiority[[1]]
[,1]       [,2]        [,3]
[1,] 0.124 0.08433333 0.031333333
[2,] 0.066 0.01200000 0.001666667
[3,] 0.810 0.90366667 0.967000000

$prob_superiority[[2]]
[,1]       [,2]
[1,] 0.1080000 0.08566667
[2,] 0.4353333 0.33966667
[3,] 0.4566667 0.57466667





### Multiple simulation example with 'Multi.BAET':
nt =3
ns=2
prior.cov = diag(25, ns*nt)
registerDoParallel(cores=3)

Multi.BAET(n.sim = 5,
           n.cores=3,
           nt=3, ns=2,
           ss.interim.es = list(c(30, 50, 100), c(20, 60, 100)),
           response.type = "gaussian",
           sig.e = 10,
           mean.response = matrix(c(seq(nt*ns)),nrow = nt, ncol = ns, byrow =F),
           prob.subpopulation = rep (1/ns, ns),
           prob.trtarm = rep (1/nt, nt),
           maxN = 100,
           upper = rep (1, ns),
           lower = rep (0, ns),
           rar = F,
           MOR = rep(0, ns),
           prob.MOR = rep(-0.5, ns),
           N.MCMC = 3000,
           prior.cov = diag(100, ns*nt),
           prior.mean = rep(0, ns*nt))

          

$est.mean
[,1]     [,2]     [,3]
[1,] 1.024811 3.225373 4.336103
[2,] 6.024678 3.610556 5.532540

$est.sd
[,1]     [,2]     [,3]
[1,] 2.302154 2.156442 2.380633
[2,] 2.255920 2.480447 2.501716

$ss.sub.dist
result.1 result.2 result.3 result.4 result.5
[1,]       52       53       51       52       55
[2,]       48       47       49       48       45

$ss.sub.mean
[1] 52.6 47.4

$ss.t.dist
result.1 result.2 result.3 result.4 result.5
100      100      100      100      100

$ss.t.mean
[1] 100

$power.sub
[1] 0 0

$computation.time
Time difference of 9.60321 secs


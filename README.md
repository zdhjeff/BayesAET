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
```
# Example:

```r
library(abind)
library(reshape2)
library(rjags)
library(doParallel)
library(MASS)
library(BayesAET)

#' @param nt number of treatment arms
#' @param ns number of subgroups
#' @param maxN the maximum sample size, trial will stop when achieving this number
#' @param nb number of subjects updated at each interim look
#' @param response.type either 'binary'(probability), 'count'(lambda) or 'gaussian'
#' @param mean.response matrix of total effect sizes: a nt (row) * ns (col) matrix, with row number indicating treatment and col number indicating subgroup
#' @param prob.subpopulation the probability of a subject coming from one specific subpopulation. default: rep (1/ns, ns)
#' @param prob.trtarm  the (initial) probability of a subject being assigned to a specific treatment arm. default: rep (1/nt, nt)
#' @param upper upper probability threshold to claim superiority (in a specific subgroup); the treatment arm wins if it goes above this threshold
#' @param lower lower probability threshold to claim superiority (in a specific subgroup); the treatment arm will be dropped if it goes below this threshold
#' @param burnin the initial sample size of before adaptation; the default is 10 times number of treatments
#' @param RAR whether using responsive adaptive randomization (RAR)
#' @param MID the minimum meaningful outcome threshold for each subgroup
#' @param prob.MID the probability threshold of being larger than the MID, treatment arms below this threshold will be dropped
#' @param N Number of MCMC samples
#' @param prior.cov the prior covariance for a multivariate normal prior
#' @param prior.mean the prior mean for a multivariate normal prior

## Gaussian outcome
nt=3
ns=2
set.seed(121)
            BAET.sim (nt=3, ns=2,
                     maxN = 500,
                     ss.interim = c(50,100,150,200,250,300,350,400,450),
                     response.type = 'gaussian',
                     sig.e = 10,
                     mean.response = matrix( c(8,12,14, 6,8,10 ), nrow = nt, ncol = ns, byrow = F),  
                     prob.subpopulation = c(0.6,0.4), # which population the patient belongs
                     prob.trtarm = rep (1/3, 3),
                     upper = rep (0.9, 2), lower = rep (0.0, 2),
                     rar = F,
                     MID = rep(6, ns),
                     prob.MID = rep(0.6, ns),
                     N = 5000,
                     prior.cov = diag(25, 3*2), prior.mean = rep(0, 3*2)
                     )   

#' @return n.interim: number of interim looks conducted to end the whole trial
#' @return trt_sub: simulated treatment arm allocation (1st column) and subgroup (2nd column)
#' @return est: treatment effect estimation(posterior mean) and the 95% Credible interval bounds for each treatment arm in each subgroup
#' @return powerind: power indicator of whether the best treatment arms is correctly selected in each subgroup
#' @return y: simlulated outcome
#' @return N_terminate: the total sample size consumed when trial ends
#' @return prob_sup_minioutcome: the probability of a treatment large than the minioutcome
#' @return prob_superiority: the probability of a treatment being the best among each subgroup

## outputs:

[1] "subgroup 1 stopped for having identified the best treatment"
[1] "subgroup 2 stopped for having identified the best treatment"

$n.interim
[1] 5

$trt_sub ## first 8 subjects
       trtarm_ind_b subpop_ind_b
  [1,]            3            1
  [2,]            3            2
  [3,]            1            1
  [4,]            3            1
  [5,]            3            1
  [6,]            3            2
  [7,]            1            1
  [8,]            2            1

$est
$est[[1]]
             trt.est  lowbound  upbound   sd.est
coef_all[1]  9.04469  6.938137 11.12762 1.077385
coef_all[2] 11.59883  9.273831 13.94035 1.169261
coef_all[3] 13.59672 11.619224 15.48728 1.001275
$est[[2]]
              trt.est  lowbound   upbound   sd.est
coef_all[4] -1.485817 -8.000615  5.235671 3.412195
coef_all[5]  3.290585 -2.829907  9.580506 3.161698
coef_all[6]  8.374851  4.739800 11.623587 1.772523

$powerind
[1] 1 1

$y ## first 8 subjects
[1]   6.86801127  14.39916139   5.94572116   8.02465238 -15.15689113   3.42099250  -1.53051597  21.55037733

$N_terminate
[1] 250

$prob_sup_minioutcome
$prob_sup_minioutcome[[1]]
          [,1]      [,2]      [,3]      [,4]      [,5]
[1,] 0.6136667 0.9603333 0.9946667 0.9953333 0.9993333
[2,] 0.8650000 0.9103333 0.9986667 1.0000000 1.0000000
[3,] 0.6616667 0.9380000 1.0000000 1.0000000 1.0000000
$prob_sup_minioutcome[[2]]
           [,1]      [,2]  [,3]      [,4]      [,5]
[1,] 0.01966667 0.0000000 0.000 0.0000000 0.0000000
[2,] 0.24066667 0.0000000 0.000 0.0000000 0.0000000
[3,] 0.70533333 0.9273333 0.907 0.9186667 0.9086667

$prob_superiority
$prob_superiority[[1]]
          [,1]      [,2]       [,3]        [,4]         [,5]
[1,] 0.2363333 0.3866667 0.07633333 0.009333333 0.0006666667
[2,] 0.5110000 0.2546667 0.19833333 0.163333333 0.0916666667
[3,] 0.2526667 0.3586667 0.72533333 0.827333333 0.9076666667
$prob_superiority[[2]]
           [,1] [,2] [,3] [,4] [,5]
[1,] 0.01966667    0    0    0    0
[2,] 0.17233333    0    0    0    0
[3,] 0.80800000    1    1    1    1





          Multi.BAET(n.sim = 10000,
                     n.cores=48,
                     nt=3, ns=2,
                     ss.interim = c(100,180,300),
                     response.type = 'gaussian',
                     sig.e = 10,
                     mean.response = matrix( c(8,12,14, 5,8,11 ), nrow = nt, ncol = ns, byrow = F),  
                     prob.subpopulation = c(0.6,0.4), # which population the patient belongs
                     prob.trtarm = rep (1/nt, nt),
                     maxN = 390,
                     upper = rep (0.9, ns), lower = rep (0.0, ns),
                     rar = F,
                     MID = rep(6, ns),
                     prob.MID = rep(0.5, ns),
                     N.MCMC = 5000,
                     prior.cov = diag(25, ns*nt), prior.mean = c(8, 12, 14, 5, 8, 11) )


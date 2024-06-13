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
                     maxN = 300,
                     nb = c(100,120,80),
                     response.type = 'gaussian',
                     sig.e = 10,
                     mean.response = matrix( c(8,12,14, 6,8,10 ), nrow = nt, ncol = ns, byrow = F),  
                     prob.subpopulation = rep (1/ns, ns), # which population the patient belongs
                     prob.trtarm = rep (1/nt, nt),
                     upper = rep (0.985, ns), lower = rep (0.0, ns),
                     rar = F,
                     MID = rep(0, ns),
                     prob.MID = rep(0.6, ns),
                     N = 3000,
                     prior.cov = diag(100, ns*nt), prior.mean = rep(0, ns*nt)
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

[1] "subgroup 1 stopped for having indentified the best treatment"
[1] "trial stopped for running out of samples"
$n.interim
[1] 3

$trt_sub ## first 10 subjects
       trtarm_ind_b subpop_ind_b
  [1,]            1            2
  [2,]            2            1
  [3,]            2            1
  [4,]            2            1
  [5,]            1            1
  [6,]            3            2
  [7,]            3            2
  [8,]            3            1
  [9,]            3            2
 [10,]            3            2

$est
$est[[1]]
              trt.est  lowbound   upbound   sd.est
coef_all[1]  6.773359  4.354134  9.230089 1.286652
coef_all[2] 10.742479  8.249738 13.317220 1.257008
coef_all[3] 15.442758 12.376759 18.755942 1.646579

$est[[2]]
              trt.est lowbound  upbound   sd.est
coef_all[4]  9.021392 6.198901 11.95611 1.435701
coef_all[5]  8.459410 5.694941 11.54432 1.498696
coef_all[6] 10.114533 7.061404 12.99339 1.490317


$powerind
[1] 1 0

$y ## first 12 subjects
  [1]  -0.09567406  15.65333353  -6.35070873  22.09899628   9.27995683   0.82785079  15.75841200  10.42735037   4.18713324  15.25805601  29.91523526   8.39749245

$N_terminate
[1] 300

$prob_sup_minioutcome
$prob_sup_minioutcome[[1]]
      [,1] [,2] [,3]
[1,] 0.995    1    1
[2,] 1.000    1    1
[3,] 1.000    1    1

$prob_sup_minioutcome[[2]]
          [,1] [,2] [,3]
[1,] 1.0000000    1    1
[2,] 1.0000000    1    1
[3,] 0.9996667    1    1


$prob_superiority
$prob_superiority[[1]]
           [,1]         [,2]  [,3]
[1,] 0.05433333 0.0003333333 0.000
[2,] 0.20633333 0.0570000000 0.007
[3,] 0.73933333 0.9426666667 0.993

$prob_superiority[[2]]
          [,1]      [,2]      [,3]
[1,] 0.4726667 0.1423333 0.2353333
[2,] 0.3336667 0.0480000 0.1476667
[3,] 0.1936667 0.8096667 0.6170000

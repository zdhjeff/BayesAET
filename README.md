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
#' @param nb number of subjects updated at each interim look
#' @param response.type either 'binary'(probability), 'count'(lambda) or 'gaussian'
#' @param totaleffect matrix of total effect sizes: a nt (row) * ns (col) matrix, with row number indicates treatment and col number indicates subgroup
#' @param prob.subpopulation the probability of a subject coming from one specific subpopulation. default: rep (1/ns, ns)
#' @param prob.trtarm  the (initial) probability of a subject being assigned to a specific treatment arm. default: rep (1/nt, nt)
#' @param maxN the maximum sample size, trial will stop when achieving this number
#' @param upper upper probability threshold to claim superiority (in a specific subgroup); the treatment arm wins if it goes above this threshold
#' @param lower lower probability threshold to claim superiority (in a specific subgroup); the treatment arm will be dropped if it goes below this threshold
#' @param burnin the initial sample size of before adaptation; the default is 10 times number of treatments
#' @param RAR whether using responsive adaptive randomization (RAR)
#' @param minioutcome the minimum meaningful outcome threshold for each subgroup
#' @param prob.minioutcome the probability threshold of being larger than the minioutcome, treatment arms below this threshold will be dropped
#' @param N Number of MCMC samples
#' @param prior.cov the prior covariance for a multivariate normal prior
#' @param prior.mean the prior mean for a multivariate normal prior

## Gaussian outcome
nt=3
ns=2

BayesAET::BAE.sim (nt=3, ns=2, nb = c(40,30,30),
                   response.type = 'gaussian',
                   totaleffect = matrix( c(1,2,2.6, 3,4,4.1 ), nrow = nt, ncol = ns, byrow = F),  ## input nt treatments for one subpopulation, then nt treatments for the second population
                   prob.subpopulation = rep (1/ns, ns), # which population the patient belongs
                   prob.trtarm = rep (1/nt, nt),
                   maxN = 100,
                   upper = rep (0.985, ns), lower = rep (0.0, ns),
                   RAR = F,
                   minioutcome = rep(3, ns),
                   prob.minioutcome = rep(0.6, ns),
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

$n.interim
[1] 3

$trt_sub # showing first 10 rows
       trtarm_ind_b subpop_ind_b
  [1,]            3            2
  [2,]            2            2
  [3,]            1            1
  [4,]            1            1
  [5,]            3            2
  [6,]            2            1
  [7,]            2            1
  [8,]            2            1
  [9,]            2            1
 [10,]            1            2


$est
$est[[1]]
             trt.est  lowbound  upbound
coef_all[1] 1.061081 0.3432057 1.823764
coef_all[2] 2.151141 1.3548472 2.942625
coef_all[3] 2.869759 2.5029568 3.229270

$est[[2]]
             trt.est lowbound  upbound
coef_all[4] 2.741041 1.793971 3.561921
coef_all[5] 3.492276 3.114791 3.895053
coef_all[6] 4.343662 4.053077 4.651546


$powerind
[1] 1 1

$y
  [1] 4.8685329 3.8876538 1.8811077 1.3981059 3.4879736 2.3411197 0.8706369 3.4330237 3.9803999 2.6327785 2.9558654 4.5697196 0.8649454 6.5016178
 [15] 2.5607600 4.6897394 1.0280022 0.2567268 3.1887923 0.1950414 4.0655549 4.1532533 4.7726117 3.4755095 3.2900536 4.7107264 2.0659024 2.7463666
 [29] 2.2914462 3.6567081 4.1011054 2.6743413 3.4104791 2.4313313 0.8648214 4.2254473 3.4993034 5.7782972 3.6874801 1.6277132 4.0253829 2.6274753
 [43] 0.9198173 3.6537509 2.9804009 4.4356172 3.0947958 2.7380527 2.4812080 2.7976843 1.5313073 3.1967868 2.8862349 4.1800917 4.0978188 2.8626455
 [57] 2.8670988 2.5962765 4.1116723 2.1243017 3.3979164 1.6259974 4.7893727 3.0441609 1.3682929 3.0431081 1.7302171 3.1893193 3.3412763 4.0685115
 [71] 4.6825860 4.0874707 3.7251452 4.4178857 3.6111944 6.7586580 5.7802782 4.8795840 4.8132405 3.4571181 4.9857784 3.6514053 2.9919454 5.9831825
 [85] 3.1710289 3.8058035 3.3850497 3.0529242 4.6989751 2.5763851 3.8938110 3.5257046 2.6098340 4.0295826 3.5691205 3.5077746 5.0811162 4.6324094
 [99] 4.0095439 4.1564905

$N_terminate
[1] 100

$prob_sup_minioutcome
$prob_sup_minioutcome[[1]]
           [,1]      [,2] [,3]
[1,] 0.00000000 0.0000000    0
[2,] 0.02366667 0.0000000    0
[3,] 0.88100000 0.1953333    0

$prob_sup_minioutcome[[2]]
          [,1]      [,2]      [,3]
[1,] 0.3140000 0.0000000 0.0000000
[2,] 0.9883333 0.9916667 0.9956667
[3,] 0.9996667 1.0000000 1.0000000


$prob_superiority
$prob_superiority[[1]]
      [,1] [,2] [,3]
[1,] 0.000    0    0
[2,] 0.006    0    0
[3,] 0.994    1    0

$prob_superiority[[2]]
            [,1]  [,2] [,3]
[1,] 0.002333333 0.000    0
[2,] 0.064666667 0.057    0
[3,] 0.933000000 0.943    1

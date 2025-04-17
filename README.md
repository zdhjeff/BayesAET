# Bayesian-Adaptive-Enrichment-Design
An R package for Bayesian adaptive enrichment trials with multiple treatment arms and multiple (categorical) subpopulations from the paper "Bayesian adaptive enrichment design in multi-arm clinical trials: The BayesAET package for R users"  <br>
Run the following command to install the package and other packages:

```r
#install the following packages if needed

library(abind)
library(reshape2)
library(rjags)
library(doParallel)
library(MASS)
library(matlib)
library(fastDummies)
library(devtools)

install_github("zdhjeff/BayesAET")

```
# Example:

```r
library(abind)
library(matlib)
library(reshape2)
library(rjags)
library(doParallel)
library(MASS)
library(fastDummies)
library(BayesAET)

#' BAET simulator: this function simulates a Bayesian adaptive enrichment trial
#' @param nt number of treatment arms
#' @param ns number of subgroups
#' @param ss.interim.es A list containing the g vectors, with gth vector indicating the accumulated number of sample sizes for each interim analysis for gth subpopulation
#' @param response.type either 'binary'(probability), 'count'(lambda) or 'gaussian'
#' @param mean.response vector of mean responses: a nt (row) * ns (col) matrix, with row number indicates treatment and col number indicates subgroup
#' @param sig.e The standard deviation of the error term when generating a gaussian outcome. The default value is 1. This parameter is only relevant when the “response.type” is set as “gaussian”.
#' @param shape The shape parameter for the gamma distribution for the standard deviation of the error term, i.e., “sig.e”. The default value is 0.001. This parameter is only relevant when the “response.type” is set as “gaussian”.
#' @param rate The rate parameter for the gamma distribution for the standard deviation of the error term, i.e., “sig.e”. The default value is 0.001. This parameter is only relevant when the “response.type” is set as “gaussian”.
#' @param prob.subpopulation the probability of a subject coming from one specific subpopulation. default: rep (1/ns, ns)
#' @param prob.trtarm  the (initial) probability of a subject being assigned to a specific treatment arm. default: rep (1/nt, nt)
#' @param maxN the maximum sample size, trial will stop when achieving this number
#' @param upper upper probability threshold to claim superiority (in a specific subgroup); the treatment arm wins if it goes above this threshold
#' @param lower lower probability threshold to claim superiority (in a specific subgroup); the treatment arm will be dropped if it goes below this threshold
#' @param rar whether using responsive adaptive randomization (rar)
#' @param rarmin.p the minimum randomization probability under rar  default:0.1
#' @param rarmax.p the maximum randomization probability under rar  default:0.9
#' @param MOR the minimum meaningful outcome threshold for each subgroup
#' @param prob.MOR the probability threshold of being larger than the MOR, treatment arms below this threshold will be dropped
#' @param N.MCMC Number of MCMC samples
#' @param prior.cov the prior covariance for a multivariate normal prior
#' @param prior.mean the prior mean for a multivariate normal prior

#' @return n.analysis: number of analysis conducted to end the whole trial
#' @return interim.sub: A vector indicating the sequence of subpopulations that reach the specified sample size threshold and trigger interim analyses.
#' @return n_terminate: the total sample size consumed when trial ends.
#' @return trt_sub: simulated treatment arm allocation (1st column) and subgroup (2nd column)
#' @return ss.sub: a vector of length ‘ns’ indicating the sample size consumed in each subpopulation.
#' @return ss.sub.trt: a ‘ns’ * ‘nt’ matrix storing the sample size for each subpopulation in each treatment arm, with the row number indicating the subpopulation and the column number indicating the treatment arm.
#' @return est: treatment effect estimation(posterior mean) and the 95% Credible interval bounds for each treatment arm in each subgroup
#' @return powerind: power indicator of whether the best treatment arms is correctly selected in each subgroup
#' @return y: simlulated outcome
#' @return prob_sup_minioutcome: the probability of a treatment large than the MOR
#' @return prob_superiority: the probability of a treatment being the best among each subgroup

## Gaussian outcome
set.seed(123)
nt =3
ns =2

BAET.sim (nt, ns,
          ss.interim.es = list(c(30, 50, 100), c(20, 60, 100)),
          response.type = "gaussian",
          sig.e = 10,
          mean.response = matrix(c(seq(nt*ns)), nrow = nt, ncol = ns, byrow =F),
          prob.subpopulation = rep (1/ns, ns),
          prob.trtarm = rep (1/nt, nt),
          maxN = 300,
          upper = rep (0.90, ns),
          lower = rep (0.10, ns),
          rar = F,
          rarmin.p = 0.1,
          rarmax.p = 0.9,
          MOR = rep(-Inf, ns),
          prob.MOR = rep(0.10, ns),
          N.MCMC = 3000,
          prior.cov = diag(100, ns*nt),
          prior.mean = rep(0, ns*nt) #notice for binary/count outcome,the prior 
          #is on log-odds/log-rate, ie, prior is put on the coef of the model
)
 

## outputs:

$n.analysis
[1] 5

$interim.sub
[1] 2 1 2 2

$ss.sub.trt

T_1 T_2 T_3
S_1   9   9  12
S_2   8 132 130

$trt_sub # only shows the first 5 data points
trtarm_ind_b subpop_ind_b
[1,]            1            2
[2,]            3            1
[3,]            1            2
[4,]            3            1
[5,]            2            1


$est
$est[[1]]
trt.est  lowbound   upbound   sd.est
coef_all[1] -0.8251961 -6.508805  5.095899 3.004195
coef_all[2]  0.1225907 -6.415069  6.725054 3.328428
coef_all[3]  8.3612976  2.910649 13.625250 2.712634

$est[[2]]
trt.est  lowbound  upbound    sd.est
coef_all[4] 0.9771983 -5.585785 7.547278 3.3739826
coef_all[5] 5.2483251  3.636972 6.908097 0.8459380
coef_all[6] 6.6165099  4.893104 8.330699 0.8727077


$powerind
[1] 1 0

$y # only shows the first 5 data points
[1]   6.15941569   6.79639483  -1.02323453  -0.33207384  -8.18575383  

$n_terminate
[1] 300

$ss.sub
[,1]
[1,]   30
[2,]  270

$prob_sup_minioutcome
$prob_sup_minioutcome[[1]]
[,1] [,2]
[1,]    1    1
[2,]    1    1
[3,]    1    1

$prob_sup_minioutcome[[2]]
[,1] [,2] [,3] [,4]
[1,]    1    1    1    1
[2,]    1    1    1    1
[3,]    1    1    1    1


$prob_superiority
$prob_superiority[[1]]
[,1] [,2]
[1,] 0.00900000    0
[2,] 0.03533333    0
[3,] 0.95566667    1

$prob_superiority[[2]]
[,1]      [,2]      [,3]      [,4]
[1,] 0.0260000 0.0000000 0.0000000 0.0000000
[2,] 0.3343333 0.1773333 0.4296667 0.1243333
[3,] 0.6396667 0.8226667 0.5703333 0.8756667








### Multiple simulation example with 'Multi.BAET':
#' Multi.BAET: This function simulates the properties of Bayesian adaptive enrichment trial
#' @param n.sim number of simulations
#' @param n.cores number of cores for computation, default is half of the total cores

#' @return est.mean: the posterior mean of each treatment (col) in each subpopulation (row); e.g. est[1,2] is treatment 2 in subpopulation 1
#' @return est.sd: the posterior sd of each treatment (col) in each subpopulation (row)
#' @return ss.sub.trt.mean: a ‘ns’ * ‘nt’ matrix storing the averaged sample size for each subpopulation in each treatment arm, with the row number indicating the subpopulation and the column number indicating the treatment arm.
#' @return ss.sub.dist: the sample size distribution for each subpopulation (row)
#' @return ss.sub.mean: the expected (mean) sample size for each subpopulation
#' @return ss.t.dist: the sample size distribution for the whole trial
#' @return ss.t.mean: the mean total sample size for the whole trial
#' @return power.sub: the power for each subpopulation
#' @return computation.time: the overall simulation time used
set.seed(123)
nt=3
ns=2

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
           rarmin.p = 0.1,
           rarmax.p = 0.9,
           MOR = rep(0, ns),
           prob.MOR = rep(0.1, ns),
           N.MCMC = 3000,
           prior.cov = diag(25, ns*nt),
           prior.mean = rep(0, ns*nt))

          
## outputs:

$est.mean
[,1]     [,2]     [,3]
[1,] 1.258367 2.397851 1.376474
[2,] 3.130782 3.975079 5.762080

$est.sd
[,1]     [,2]     [,3]
[1,] 2.400057 2.171975 2.074618
[2,] 2.209861 2.194168 2.094379

$ss.sub.trt.mean

T_1  T_2  T_3
S_1 13.6 17.8 18.0
S_2 15.8 16.2 18.6

$ss.sub.dist
result.1 result.2 result.3 result.4 result.5
[1,]       47       52       48       50       50
[2,]       53       48       52       50       50

$ss.sub.mean
[1] 49.4 50.6

$ss.t.dist
result.1 result.2 result.3 result.4 result.5 
100      100      100      100      100 

$ss.t.mean
[1] 100

$power.sub
[1] 0 0

$computation.time
Time difference of 7.050495 secs





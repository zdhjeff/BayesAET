# Bayesian-Adaptive-Enrichment-Design
An R package for Bayesian adaptive enrichment trials with multiple treatment arms and multiple pre-specified (categorical) subpopulations from the paper: <br>
"Zhan D, Ouyang Y, Vila-Rodriguez F, Karim ME, Wong H. Bayesian adaptive enrichment design in multi-arm clinical trials: The BayesAET package for R users. Computer Methods and Programs in Biomedicine. 2025 May 5:108833." https://www.sciencedirect.com/science/article/pii/S0169260725002500 <br>
Run the following command to install the package and other packages:

```r
# Install the following packages if needed

library(abind)
library(reshape2)
library(rjags)
library(doParallel)
library(MASS)
library(matlib)
library(fastDummies)
library(devtools)

# install the BayesAET package
install_github("zdhjeff/BayesAET")

```
# Example:

```r
library(abind)
library(reshape2)
library(rjags)
library(doParallel)
library(MASS)
library(matlib)
library(fastDummies)
library(BayesAET)

#' BAET simulator: this function simulates a Bayesian adaptive enrichment trial
#' @param nt The number of treatment arms.
#' @param ns The number of subpopulations.
#' @param ss.interim.es A list containing g vectors, with g-th vector indicating the accumulated sample sizes at each interim analysis for g-th subpopulation.
#' @param response.type The outcome type, one of “gaussian”, “binary” or “count”.
#' @param mean.response A ‘nt’ * ‘ns’  matrix of the assumed true mean responses for each treatment (row number) in each subpopulation (column number). For binary and count outcomes, the inputs are the probabilities and rates, respectively.
#' @param sig.e The standard deviation of the error term when generating a gaussian outcome. The default value is 1. This parameter is only relevant when the “response.type” is set as “gaussian”.
#' @param shape The shape parameter for the gamma distribution for the standard deviation of the error term, i.e., “sig.e”. The default value is 0.001. This parameter is only relevant when the “response.type” is set as “gaussian”.
#' @param rate The rate parameter for the gamma distribution for the standard deviation of the error term, i.e., “sig.e”. The default value is 0.001. This parameter is only relevant when the “response.type” is set as “gaussian”.
#' @param prob.subpopulation A vector of length ‘ns’ indicating the probabilities of a participant belonging to each of the subpopulations. The default is an equal probability for each subpopulation.
#' @param prob.trtarm  A vector of length ‘nt’ indicating the (initial) probability of a participant being assigned to each treatment arm. The default is an equal probability to each treatment arm.
#' @param maxN The maximum sample size.
#' @param upper A vector of length ‘ns’ indicating the posterior probability threshold above which a treatment will be declared best (superiority) for each subpopulation. The default value is 0.90.
#' @param lower A vector of length ‘ns’ indicating the posterior probability threshold below which a treatment will be dropped for each subpopulation. The default value is 0.10.
#' @param rar A logical indicator of whether using responsive adaptive randomization (RAR). The default setting is ‘false’, i.e., balanced randomization.
#' @param rarmin.p The minimum randomization probability under ‘rar’. The default value is 0.1.
#' @param rarmax.p The maximum randomization probability under ‘rar’. The default value is 0.9.
#' @param MOR A vector of length ‘ns’ indicating the minimum outcome requirement threshold for each subpopulation. The default value is ‘-Inf’.
#' @param prob.MOR A vector of length ‘ns’ indicating the posterior probability threshold below which a treatment will be dropped due to the outcome not achieving the ‘MOR’ for each subpopulation. The default value is 0.10.
#' @param N.MCMC The number of MCMC samples (excluding burn in samples)
#' @param prior.cov A (‘nt’ * ‘ns’) by (‘nt’ * ‘ns’) square matrix for the prior covariance of a multivariate normal prior for the mean responses. The default is a diagonal matrix with ‘100’ on the diagonal.
#' @param prior.mean A vector of length ‘nt’ * ‘ns’ for prior mean of a multivariate normal prior for the mean responses. The default is a vector of zeros.


#' @return n.analysis: The total number of analyses, i.e., the number of interim analysis plus the final analysis.
#' @return interim.sub: A vector indicating the sequence of subpopulations that reach the specified sample size threshold and trigger interim analyses.
#' @return n_terminate: The total sample size consumed when the trial ended.
#' @return trt_sub: A two-column matrix with the first column displaying the treatment each participant received and the second column displaying the subpopulation in which the participant belongs.
#' @return ss.sub: A vector of length ‘ns’ indicating the sample size consumed in each subpopulation.
#' @return ss.sub.trt: A ‘ns’ * ‘nt’ matrix storing the sample size for each subpopulation in each treatment arm, with the row number indicating the subpopulation and the column number indicating the treatment arm.
#' @return est: A list of length ‘ns’ with each component summarizing the posterior distribution for the mean responses within a subpopulation at trial termination. Each component is a ‘nt’ * 4 data frame with rows representing treatments and columns displaying the posterior mean, the lower and upper bound of 95% Credible Interval, and the SD of the posterior distribution, respectively.
#' @return powerind: A vector of length ‘ns’ indicating whether the best treatment arm was correctly selected in each subpopulation at trial termination. The value is ‘1’ if correctly selected and ‘0’ otherwise.
#' @return y: A vector for the simulated outcomes.
#' @return prob_sup_minioutcome: A list of length ‘ns’ with each component summarizing the probability of superiority over the fixed MOR within a subpopulation. Each component is a matrix with rows representing treatments and columns displaying the probability of superiority over the fixed MOR at each interim look.
#' @return prob_superiority: A list of length ‘ns’ with each component summarizing the probability of superiority among different treatments within a subpopulation. Each component is a matrix with rows indicating treatments and columns displaying the probability of superiority for each treatment at each interim look.


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
#' @param n.sim The number of simulations
#' @param n.cores The number of cores used for computation, default is half of the total cores

#' @return est.mean: A ‘ns’ * ‘nt’ matrix storing the mean of the posterior mean, with row number indicating subpopulation and column number indicating treatment arm
#' @return est.sd: A ‘ns’ * ‘nt’ matrix storing the mean of the posterior SD, with row number indicating subpopulation and column number indicating treatment arm
#' @return ss.sub.trt.mean: A ‘ns’ * ‘nt’ matrix storing the averaged sample size for each subpopulation in each treatment arm, with the row number indicating the subpopulation and the column number indicating the treatment arm.
#' @return ss.sub.dist: A ‘ns’ * ‘n.sim’ matrix storing the sample size consumed for each subpopulation in each run, with row number indicating subpopulation and column number indicating the specific run
#' @return ss.sub.mean: A vector length ‘ns’ indicating the average sample size consumed in each subpopulation
#' @return ss.t.dist: A ‘n.sim’ vector storing the overall sample size consumed in each run
#' @return ss.t.mean: A number of average overall sample size consumed
#' @return power.sub: A vector length ‘ns’ indicating the (empirical) power estimated in each subpopulation
#' @return computation.time: Time consumed for computation

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





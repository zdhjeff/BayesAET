# Bayesian-Adaptive-Enrichment-Design
An R package for Bayesian adaptive enrichment trials with multiple treatment arms and multiple (categorical) subpopulations  <br>
Run the following command to install the package:

```r
library(devtools) 
install_github("zdhjeff/BayesAET")
```
# Example:

```r
library(abind)
library(reshape2)
library(rjags)
library(doParallel)
library(MASS)


## Gaussian outcome
nt=3
ns=2

BayesAET::BAE.sim (nt=3, ns=2, nb = 35,
                   response.type = 'gaussian',
                   totaleffect = matrix( c(1,2,2.6, 3,4,4.1 ), nrow = nt, ncol = ns, byrow = F),  ## input nt treatments for one subpopulation, then nt treatments for the second population
                   prob.subpopulation = rep (1/ns, ns), # which population the patient belongs
                   prob.trtarm = rep (1/nt, nt),
                   maxN = 100,
                   upper = rep (0.985, ns), lower = rep (0.0, ns),
                   burnin = 10*nt,
                   RAR = F,
                   minioutcome = rep(3, ns),
                   prob.minioutcome = rep(0.6, ns),
                   N = 3000,
                   prior.cov = diag(100, ns*nt), prior.mean = rep(0, ns*nt)
                   )

## outputs:

$n.interim
[1] 3

$trt_sub # showing first 10 rows
       trtarm_ind_b subpop_ind_b
  [1,]            2            2
  [2,]            2            1
  [3,]            3            1
  [4,]            3            2
  [5,]            3            2
  [6,]            1            2
  [7,]            2            1
  [8,]            2            1
  [9,]            3            2
 [10,]            1            2

$est
$est[[1]]
             trt.est  lowbound  upbound
coef_all[1] 1.189086 0.3632326 1.985910
coef_all[2] 2.456063 1.8147152 3.059142
coef_all[3] 2.465858 1.6283231 3.273285

$est[[2]]
             trt.est lowbound  upbound
coef_all[4] 2.903149 2.229916 3.568666
coef_all[5] 3.858213 3.543193 4.214175
coef_all[6] 4.150933 3.871708 4.422605


$powerind
[1] 0 0

$y
  [1] 4.2516860 3.2920629 2.3805870 3.5796725 4.6582085 2.6744351 2.4543218 3.0837509 4.1367432 1.9045497 4.4182898 0.2001504 1.3880218 0.5733419
 [15] 2.6797310 3.1844862 4.1864996 3.0444913 3.7310783 2.2125619 3.5149300 2.7700670 3.0137887 1.1878443 4.8039606 2.3466388 0.9688889 1.6248045
 [29] 2.2266699 3.7278098 4.0998005 2.4607859 1.4060860 1.4670152 4.1498216 4.8530366 4.1281762 5.0928433 4.6904302 1.6578910 2.6550487 3.6910025
 [43] 2.7251260 3.4352137 6.1333514 3.8557043 4.1013414 4.6495323 4.0643813 4.3996691 5.2897185 4.5231373 3.8403787 4.5751224 3.5342806 3.6607689
 [57] 2.9966185 6.2590030 1.9143114 3.3288321 1.5792420 4.6572712 4.9989605 5.0332766 2.5669356 4.7395060 3.3000157 3.1156482 2.6175842 3.3006880
 [71] 4.2638275 4.1287950 4.7658824 4.5517657 3.9631188 4.3070268 4.1254080 2.5950106 5.9697374 3.4289164 4.5274576 4.9203605 4.9414372 5.2598798
 [85] 3.6077723 3.9808137 3.4195916 4.0210006 3.7728704 3.7622788 4.5739454 3.2587985 4.2441389 3.8907285 3.8396958 4.3269594 5.1560482 4.3563673
 [99] 4.0202971 4.5702877

$N_terminate
[1] 100

$prob_sup_minioutcome
$prob_sup_minioutcome[[1]]
           [,1] [,2] [,3]
[1,] 0.00000000    0    0
[2,] 0.04300000    0    0
[3,] 0.04433333    0    0

$prob_sup_minioutcome[[2]]
          [,1]      [,2] [,3]
[1,] 0.3276667 0.0000000    0
[2,] 0.7233333 0.9986667    1
[3,] 1.0000000 1.0000000    1


$prob_superiority
$prob_superiority[[1]]
          [,1] [,2] [,3]
[1,] 0.0000000    0    0
[2,] 0.5216667    0    0
[3,] 0.4783333    0    0

$prob_superiority[[2]]
            [,1]  [,2]  [,3]
[1,] 0.001333333 0.000 0.000
[2,] 0.043000000 0.065 0.092
[3,] 0.955666667 0.935 0.908


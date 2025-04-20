# load the required packages
library(ggplot2)
library(abind)
library(reshape2)
library(rjags)
library(doParallel)
library(MASS)
library(roxygen2)
library(fastDummies)
library(parallel)

#' @import abind
#' @import reshape2
#' @import rjags
#' @import doParallel
#' @import foreach
#' @import MASS
#' @import roxygen2
#' @import fastDummies


# Helper function to adjust probability under rar:
adjust_prob <- function(prob, rarmin.p, rarmax.p) {
  # Step 1: Handle special cases
  if (any(is.na(prob))) stop("Input probability vector contains NA values")
  #if (sum(prob) == 0) return(rep(1 / length(prob), length(prob)))  # Assign equal probabilities

  # Step 2: Normalize the probability vector
  prob <- prob / sum(prob)

  # Step 3: Identify values that exceed limits
  over_max <- prob > rarmax.p
  under_min <- (prob < rarmin.p & prob > 0)

  # Step 4: Fix probabilities that exceed limits
  prob[over_max] <- rarmax.p
  prob[under_min] <- rarmin.p

  # Step 5: Redistribute remaining probability to maintain sum = 1
  remaining_prob <- 1 - sum(prob)  # The total amount we need to redistribute
  free_indices <- !(over_max | under_min)  # Indices that can be adjusted

  if (sum(free_indices) > 0 && sum(prob[free_indices]) > 0) {
    # Redistribute remaining probability proportionally
    prob[free_indices] <- prob[free_indices] + remaining_prob * (prob[free_indices] / sum(prob[free_indices]))
  } else if (sum(over_max) > 0) {
    # If no free elements, force sum to 1 by slightly adjusting max values
    prob[over_max] <- prob[over_max] + (remaining_prob / sum(over_max))
  }

  return(prob)
}


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

#' @examples
#' nt =3
#' ns =2
#' BAET.sim (nt=3, ns=2,
#'          ss.interim.es = list(c(30, 50, 100), c(20, 60, 100)),
#'          response.type = "gaussian",
#'          sig.e = 10,
#'          mean.response = matrix(c(seq(nt*ns)), nrow = nt, ncol = ns, byrow =F),
#'          prob.subpopulation = rep (1/ns, ns),
#'          prob.trtarm = rep (1/nt, nt),
#'         maxN = 300,
#'          upper = rep (0.90, ns),
#'          lower = rep (0.10, ns),
#'         rar = F,
#'          rarmin.p = 0.1,
#'          rarmax.p = 0.9,
#'          MOR = rep(-Inf, ns),
#'          prob.MOR = rep(0.10, ns),
#'          N.MCMC = 3000,
#'          prior.cov = diag(100, ns*nt),
#'          prior.mean = rep(0, ns*nt) #notice for binary/count outcome,the prior
#'          #is on log-odds/log-rate, ie, prior is put on the coef of the model
#')

#' @export


BAET.sim = function(nt, ns,
                    ss.interim.es,
                    response.type = "gaussian",
                    sig.e = 10,
                    mean.response = matrix(c(seq(nt*ns) ), nrow = nt, ncol = ns, byrow = F),  ## input nt treatments for one subpopulation, then nt treatments for the second population
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
                    prior.mean = rep(0, ns*nt) ## notice for binary/count outcome,the prior is on log-odds/log-rate, ie, prior is put one the coef of the model
){

  # now begin the functions:
  response.types = c("gaussian", "binary", "count")  # set allowed outcome types
  if (!(response.type %in% response.types)) return(
    "Allowed outcome distributions are \"gaussian\", \"poisson\", or \"binomial\".")

  supority_comp = function(x) {
    n = length(x)
    check = c()
    for (i in 1:n) check[i] = (x[i] == max(x))
    return(check)
  }


  y = NULL # outcome
  x_trtarm = NULL
  x_subpop = NULL
  trtarm_ind = NULL
  subpop_ind = NULL
  trt_sub = NULL
  x_design = NULL

  # Initialize counts for each category
  counts <- rep(0, ns)
  # Create a list to hold all vectors
  threshold_list <- ss.interim.es
  # Initialize threshold indices for each group (starting at the first threshold)
  threshold_indices <- rep(1, ns)
  # Set the initial thresholds for each group (first value in each vector)
  current_thresholds <- sapply(threshold_list, function(v) v[1])
  # Storage for generated data
  generated_data <- NULL
  total_samples <- 0  # Counter for total samples generated


  check = list() # check results for each subpopualtion (ns), for the nt treatments
  thetaall<-list()
  prob_superiority = list ()# posterior probability outputs, in each ns subpopulation
  prob_assign = list()
  prob.subpop=prob.subpopulation
  treat_prob = list()
  treat_result = list()

  for (i in 1:ns) {
    check[[i]] = array(0, dim = c(N.MCMC, nt)) # check the which MCMC samples are the largest, and return it to 1 for that treatment
    thetaall[[i]] = array (rep(0,N.MCMC*nt), dim = c(N.MCMC, nt))
    prob_superiority[[i]] = array(rep(1/nt, nt), dim = c(nt, 1))
    prob_assign[[i]] =  prob.trtarm
    treat_prob[[i]] = rep(1,nt)
    treat_result[[i]] = rep(1,nt)
  }
  treat_result_record = treat_result # the records for each interim look of treat_result
  prob_sup_minioutcome = treat_prob ## the records for each interim look of treatprob

  ntj = rep(0,ns)
  powerind = rep(0, ns) # power indicator for each subpopulation
  C1 = rep(0, ns) # Creteria 1 indicator for each subpopulation
  C2 = rep(0, ns) # Creteria 2 indicator for each subpopulation
  earlystop = rep(0,ns)

  trt.est = list()
  lowbound=list()
  upbound = list()
  sd.est = list()
  est= list()

  nb= c()
  j = 0
  sg = c() # indicator of which subgroup is checked during interim looks
  #######################  Start the loop ###########################
  repeat{
    j = j + 1
    while (any(threshold_indices <= sapply(threshold_list, length)) && total_samples < maxN) {

      # Generate a single observation from the multinomial distribution
      new_obs <- rmultinom(1, size = 1, prob = prob.subpop)

      # Update counts for each category
      counts <- counts + new_obs
      total_samples <- total_samples + 1

      # Store the generated observation
      generated_data <- rbind(generated_data, t(new_obs))

      if (any(counts >= current_thresholds)) {
        break
      }

    }
    # Check if any of the counts have reached their current threshold
    for (i in 1:ns) {
      if (counts[i] >= current_thresholds[i] && threshold_indices[i] <= length(threshold_list[[i]])) {
        sg[j] = i
        # Move to the next threshold in the vector for that category if it's not the last threshold
        if (threshold_indices[i] < length(threshold_list[[i]])) {
          threshold_indices[i] <- threshold_indices[i] + 1
          current_thresholds[i] <- threshold_list[[i]][threshold_indices[i]]
        } else {
          # Set threshold to a very large number if the last threshold is reached
          current_thresholds[i] <- Inf
        }
      }
    }


    # if (sum(nb[1:j]) > maxN) nb[j] = maxN - sum(nb[1:(j-1)])

    x_subpop_b = generated_data # this is the matrix indicating which subpopulation for the patient, in the batch
    nb[j] = nrow(x_subpop_b) # this is the batch size
    generated_data =  NULL # reset this matrix after giving its value to x_subpop_b
    x_trtarm_b = matrix( rep(0,nb[j]*nt), nrow = nb[j], ncol=nt)

    for (i in 1:ns) {
      if (sum(prob_assign[[i]])!= 0){
        x_trtarm_b[which(x_subpop_b[,i]==1) ,] = t(rmultinom(length(which(x_subpop_b[,i]==1)), 1, prob = prob_assign[[i]]))
      }

    }# this is the matrix indicating which treatment arm for the patient, in the batch


    trtarm_ind_b = apply( x_trtarm_b, 1, function(z) which(z==1)) # this is the treatment indicatior for each nb patient, in the batch
    subpop_ind_b = apply( x_subpop_b, 1, function(z) which(z==1)) # this is the subpopulation indicatior for each nb patient, in the batch
    trt_sub_b = cbind(trtarm_ind_b, subpop_ind_b) # this is the nb* 2 matrix with first col shows the treatment arm, 2nd col shows the subpop, in the batch

    design_b = matrix(rep(0, nt*ns*nb[j]), nrow = nb[j], ncol = (nt*ns))
    for (i in 1:nb[j]) {
      indmatrix = matrix(rep(0, ns*nt), nrow = ns, ncol = nt)
      indmatrix[trt_sub_b[i,2],trt_sub_b[i,1]] =1
      design_b[i,] = c(t(indmatrix))

    }

    if (response.type == 'gaussian') {
      yb = apply (trt_sub_b , 1, function(z) mean.response[z[1], z[2]]) + rnorm(nb[j],0,sig.e)
    }
    if (response.type == 'binary') {
      yb = apply (trt_sub_b , 1, function(z) rbinom(1,1,mean.response[z[1], z[2]] )  )
    }
    if (response.type == 'count') {
      yb = apply (trt_sub_b , 1, function(z) rpois(1, mean.response[z[1], z[2]] )  )
    }


    x_trtarm = abind(x_trtarm, x_trtarm_b , along=1 )
    x_subpop = abind(x_subpop, x_subpop_b , along=1 )
    x_design = rbind(x_design, design_b)
    trtarm_ind = abind(trtarm_ind, trtarm_ind_b , along=1 )
    subpop_ind = abind(subpop_ind, subpop_ind_b , along=1 )
    trt_sub = abind(trt_sub, trt_sub_b , along=1 )
    y = c(y, yb)

    ################ now begin the JAGS script ########################
    np = nt * ns

    jags_works = 0
    attempt = 1
    if(response.type=='gaussian'){
      while (jags_works == 0){
        # Following call should work will with all outcome types; not all parameters passed
        # are needed for every outcome type

        gaussian.model<-"model{

          for ( j in 1:p ){

            y[j] ~ dnorm(mu[j], sig.e^2)

            mu[j] <- (x_design[j,] %*% coef_all)
          }

          # priors
          coef_all[1:np] ~ dmnorm.vcov(Mean, prior.cov)
          sig.e ~ dgamma(0.0001, 0.0001)

        }"


        jags=try(jags.model(textConnection(gaussian.model),
                            data = list("p" = length(y), "y" = y, "x_design" = x_design,  "np" = np,
                                        #"nt" =nt,
                                        "prior.cov"= prior.cov, "Mean" = prior.mean),
                            inits = list( # "theta_all[1:np]"= rep (1, np),
                              "sig.e" = 1),
                            quiet=T) )
        jags_works = as.numeric(length(jags) > 1 | attempt == 10)
        attempt = attempt + 1
      }
      update(jags, 5500) # more burn-in
      my.samples = coda.samples(jags, variable.names = "coef_all", n.iter = N.MCMC)
      samples = my.samples[[1]]
      thetaall_new<- samples # theta_new stores the total effect, with seq t1s1,t2s1,t3s1, t1s2,t2s2,t3s2...
      thetaall_new<-as.data.frame(thetaall_new)

    }


    if(response.type=='binary'){
      while (jags_works == 0){
        # Following call should work will with all outcome types; not all parameters passed
        # are needed for every outcome type

        binary.model<-"model{

          for ( j in 1:p ){

            y[j] ~ dbern(mu[j])

            logit(mu[j]) <- (x_design[j,] %*% coef_all)
          }

          # priors
          coef_all[1:np] ~ dmnorm.vcov(Mean, prior.cov)

        }"


        jags=try(jags.model(textConnection(binary.model),
                            data = list("p" = length(y), "y" = y, "x_design" = x_design,  "np" = np,
                                        #"nt" =nt,
                                        "prior.cov"= prior.cov, "Mean" = prior.mean),
                            # inits = list( # "theta_all[1:np]"= rep (1, np),
                            #   "sig.e" = 1),
                            quiet=T) )
        jags_works = as.numeric(length(jags) > 1 | attempt == 10)
        attempt = attempt + 1
      }
      update(jags, 5500) # more burn-in
      my.samples = coda.samples(jags, variable.names = "coef_all", n.iter = N.MCMC)
      samples = my.samples[[1]]
      thetaall_new<- samples # theta_new stores the total effect but in log-odd, with seq t1s1,t2s1,t3s1, t1s2,t2s2,t3s2...
      thetaall_new<- exp(thetaall_new)/(exp(thetaall_new) +1) # convert the log-odd to probability
      thetaall_new<-as.data.frame(thetaall_new)
    }

    if(response.type=='count'){
      while (jags_works == 0){
        # Following call should work will with all outcome types; not all parameters passed
        # are needed for every outcome type

        count.model<-"model{

          for ( j in 1:p ){

            y[j] ~ dpois(mu[j])

            log(mu[j]) <- (x_design[j,] %*% coef_all)
          }

          # priors
          coef_all[1:np] ~ dmnorm.vcov(Mean, prior.cov)

        }"


        jags=try(jags.model(textConnection(count.model),
                            data = list("p" = length(y), "y" = y, "x_design" = x_design,  "np" = np,
                                        #"nt" =nt,
                                        "prior.cov"= prior.cov, "Mean" = prior.mean),
                            # inits = list( # "theta_all[1:np]"= rep (1, np),
                            #   "sig.e" = 1),
                            quiet=T) )
        jags_works = as.numeric(length(jags) > 1 | attempt == 10)
        attempt = attempt + 1
      }
      update(jags, 5500) # more burn-in
      my.samples = coda.samples(jags, variable.names = "coef_all", n.iter = N.MCMC)
      samples = my.samples[[1]]
      thetaall_new<- samples # theta_new stores the total effect but in log-rate, with seq t1s1,t2s1,t3s1, t1s2,t2s2,t3s2...
      thetaall_new<- exp(thetaall_new) # convert the log-rate to rate
      thetaall_new<-as.data.frame(thetaall_new)
    }


    grouping <- rep(1:ns, each = nt)
    thetaall_new_split <- split.default(thetaall_new, f = grouping) # split the trts into different subpopulations
    thetaall_new_split <- lapply(thetaall_new_split, as.data.frame)

    for (i in 1:ns){
      thetaall[[i]] = abind(thetaall[[i]], thetaall_new_split[[i]], along = 3) # theta stores the total effect for each cycle, with seq t1s1,t2s1,t3s1, t1s2,t2s2,t3s2...
    } # note that the first cohort is 0, so when comparing the j interim look, use j+1 data


    ##### The above sampling process is done, the following is the decision part, based on the posterior probability ####
    ### this for loop is for decisions in each subpopulation of ns
    if (length(y) == maxN){
      for (i in 1:ns){

        if(sum(prob_assign[[i]] == 0) >0) { # this means there is at least on treatment has a probability of assignment equals to 0
          check[[i]][,which(prob_assign[[i]] == 0)] =0 # this indicates that once the probability of assigning subject to a treatment is 0, then the check for that trt is 0
        }

        if(sum(prob_assign[[i]]>0) ==1){ # this means only one treatment left
          check[[i]][,which(prob_assign[[i]]>0)] = rep(1,N.MCMC)
        }
        else if(sum(prob_assign[[i]]>0) > 1){
          check[[i]][,which(prob_assign[[i]]>0)] = t(apply( (thetaall[[i]])[,which(prob_assign[[i]]>0),j+1,drop=F], 1, supority_comp))# sup_check is function to check which one is largest and then return it to 1
        }

        prob_superiority[[i]] = abind(prob_superiority[[i]], apply(check[[i]], 2, mean), along = 2)

        treat_prob[[i]][which(treat_result[[i]]==0)]=0 ## once the subgroup failed in mini check, it is zero
        treat_prob[[i]][which(treat_result[[i]]==1)]= as.vector( (apply((thetaall_new_split[[i]][,which(treat_result[[i]]==1),drop=F]-MOR[i] >=0 ),2,sum))/N.MCMC ) ## minioutcome_1 is the minimum outcome needed, treatp is the probablity of larger than minioutcome_1
        treat_result[[i]]=as.numeric(treat_prob[[i]]>= prob.MOR[i] ) ## whether that the treatment arm is larger than the threshold
        prob_assign[[i]][which(treat_result[[i]]==0)]=0 # this means once the trt fails in mini check, that trt stops receiving subjects

        prob_sup_minioutcome[[i]]=abind(prob_sup_minioutcome[[i]],treat_prob[[i]],along = 2) ## this is only to record treat_prob[[i]]
        treat_result_record[[i]]= abind(treat_result_record[[i]],treat_result[[i]],along = 2)## this is only to record treat_result[[i]]

        ntj[i] = length(which(apply(check[[i]], 2, mean)>lower[i] & treat_result[[i]]==1)) ## No. of groups left
        if (ntj[i]!=0){
          prob_assign[[i]] = rep(0, nt)
          if (rar == T  & ntj[[i]]>1)  {
            prob_assign[[i]][which( apply(check[[i]], 2, mean)>lower[i] & treat_result[[i]]==1)] =
              ((ntj[i] - 1)/ntj[i])*sqrt(prob_superiority[[i]][which(apply(check[[i]], 2, mean)>lower[i]& treat_result[[i]]==1),threshold_indices[i]])
           # adjust the probability to make sure the probs are within a specified range
             prob_assign[[i]][which( apply(check[[i]], 2, mean)>lower[i] & treat_result[[i]]==1 & !0)] = adjust_prob(c(prob_assign[[i]][which( apply(check[[i]], 2, mean)>lower[i] & treat_result[[i]]==1)]),rarmin.p = rarmin.p, rarmax.p = rarmax.p )
          } else prob_assign[[i]][which( apply(check[[i]], 2, mean)>lower[i]& treat_result[[i]]==1)] = 1/ntj[i]

        }

        ## C1 first stopping criteria
        if (max(apply(check[[i]], 2, mean)) >= upper[i] & prob_assign[[i]][which.max(apply(check[[i]], 2, mean))] != 0 ){ # the posterior should larger than threshold and that arm is still in

          #prob_assign[[i]] = rep(0, nt)
          prob.subpop[i]=0

          C1[i]=1 ## this indicates the best subgroup is found

          if (which.max(apply(check[[i]], 2, mean)) == which.max(mean.response[1:nt,i])) powerind[i]=1 # this indicates the best subgroup was found and it is truely the best

        }

        ## C2 second stopping criteria
        if (sum(prob_assign[[i]]==0)==nt) { # this means all trts in the ith subpopulation are failed
          prob_assign[[i]] = rep(0, nt)
          prob.subpop[i]=0
          ########

          C2[i]=1 ## this indicates all subgroups failed
        }
      }
    }else{
      i = sg[j] # only make decisions for this subgroup
      if(sum(prob_assign[[i]] == 0) >0) { # this means there is at least on treatment has a probability of assignment equals to 0
        check[[i]][,which(prob_assign[[i]] == 0)] =0 # this indicates that once the probability of assigning subject to a treatment is 0, then the check for that trt is 0
      }

      if(sum(prob_assign[[i]]>0) ==1){ # this means only one treatment left
        check[[i]][,which(prob_assign[[i]]>0)] = rep(1,N.MCMC)
      }
      else if(sum(prob_assign[[i]]>0) > 1){
        check[[i]][,which(prob_assign[[i]]>0)] = t(apply( (thetaall[[i]])[,which(prob_assign[[i]]>0),j+1,drop=F], 1, supority_comp))# sup_check is function to check which one is largest and then return it to 1
      }

      prob_superiority[[i]] = abind(prob_superiority[[i]], apply(check[[i]], 2, mean), along = 2)

      treat_prob[[i]][which(treat_result[[i]]==0)]=0 ## once the subgroup failed in mini check, it is zero
      treat_prob[[i]][which(treat_result[[i]]==1)]= as.vector( (apply((thetaall_new_split[[i]][,which(treat_result[[i]]==1),drop=F]-MOR[i] >=0 ),2,sum))/N.MCMC ) ## minioutcome_1 is the minimum outcome needed, treatp is the probablity of larger than minioutcome_1
      treat_result[[i]]=as.numeric(treat_prob[[i]]>= prob.MOR[i] ) ## whether that the treatment arm is larger than the threshold
      prob_assign[[i]][which(treat_result[[i]]==0)]=0 # this means once the trt fails in mini check, that trt stops receiving subjects

      prob_sup_minioutcome[[i]]=abind(prob_sup_minioutcome[[i]],treat_prob[[i]],along = 2) ## this is only to record treat_prob[[i]]
      treat_result_record[[i]]= abind(treat_result_record[[i]],treat_result[[i]],along = 2)## this is only to record treat_result[[i]]

      ntj[i] = length(which(apply(check[[i]], 2, mean)>lower[i] & treat_result[[i]]==1)) ## No. of groups left
      if (ntj[i]!=0){
        prob_assign[[i]] = rep(0, nt)
        if (rar == T  & ntj[[i]]>1)  {
          prob_assign[[i]][which( apply(check[[i]], 2, mean)>lower[i] & treat_result[[i]]==1)] =
            ((ntj[i] - 1)/ntj[i])*sqrt(prob_superiority[[i]][which(apply(check[[i]], 2, mean)>lower[i]& treat_result[[i]]==1),threshold_indices[i]])
          # adjust the probability to make sure the probs are within a specified range
          prob_assign[[i]][which( apply(check[[i]], 2, mean)>lower[i] & treat_result[[i]]==1 & !0)] = adjust_prob(c(prob_assign[[i]][which( apply(check[[i]], 2, mean)>lower[i] & treat_result[[i]]==1)]),rarmin.p = rarmin.p, rarmax.p = rarmax.p )
        } else prob_assign[[i]][which( apply(check[[i]], 2, mean)>lower[i]& treat_result[[i]]==1)] = 1/ntj[i]

      }

      ## C1 first stopping criteria
      if (max(apply(check[[i]], 2, mean)) >= upper[i] & prob_assign[[i]][which.max(apply(check[[i]], 2, mean))] != 0 ){ # the posterior should larger than threshold and that arm is still in

        #prob_assign[[i]] = rep(0, nt)
        prob.subpop[i]=0

        C1[i]=1 ## this indicates the best subgroup is found

        if (which.max(apply(check[[i]], 2, mean)) == which.max(mean.response[1:nt,i])) powerind[i]=1 # this indicates the best subgroup was found and it is truely the best

      }

      ## C2 second stopping criteria
      if (sum(prob_assign[[i]]==0)==nt) { # this means all trts in the ith subpopulation are failed
        prob_assign[[i]] = rep(0, nt)
        prob.subpop[i]=0
        ########

        C2[i]=1 ## this indicates all subgroups failed
      }


    }


    for (i in 1:ns) {
      earlystop[i] = (C1[i]==1 | C2[i] ==1) # this means the stop sign for each subpopulation
    }

    if ( sum(earlystop)==ns | length(y) >= maxN) {

      # for (i in 1:ns) {
      #   if (C1[i]==1) print(sprintf("subgroup %d stopped for having indentified the best treatment", i))
      #   if (C2[i]==1) print(sprintf("subgroup %d stopped for futility of all treatments", i))
      # }
      # if (length(y) >= maxN) print("trial stopped for running out of samples")
      break

    }  # this means if all subpopulation meet stop criterion , or sample size is reached, stop the trial




  }
  #### now summarize the outputs ###

  for (i in 1:ns) {
    trt.est[[i]] = apply(thetaall[[i]][,,j+1],2,mean) # j+1 because the first is created to store data.
    sd.est[[i]] = apply(thetaall[[i]][,,j+1],2,sd)
    lowbound[[i]] =apply(thetaall[[i]][,,j+1],2,function(x) quantile(x,0.025))
    upbound[[i]] =apply(thetaall[[i]][,,j+1],2,function(x) quantile(x,0.975))
    est[[i]] = data.frame(trt.est = trt.est[[i]], lowbound = lowbound[[i]], upbound = upbound[[i]], sd.est = sd.est[[i]] )

    prob_sup_minioutcome[[i]] = prob_sup_minioutcome[[i]][,-1]
    treat_result_record[[i]] =treat_result_record[[i]][,-1]
    prob_superiority[[i]] = prob_superiority[[i]][,-1]
  }
  # Create a summary table for trt_sub
    ss.sub.trt <- table(trt_sub[, "subpop_ind_b"], trt_sub[, "trtarm_ind_b"])
    ss.sub.trt <- as.matrix(ss.sub.trt)
    rownames(ss.sub.trt) <- paste0("S_", rownames(ss.sub.trt))
    colnames(ss.sub.trt) <- paste0("T_", colnames(ss.sub.trt))

  out = list(
    n.analysis =j,
    interim.sub = sg,
    ss.sub.trt =ss.sub.trt,
    trt_sub =trt_sub,
    est =est,
    powerind = powerind,
    y =y,
    n_terminate = length(y),
    ss.sub=counts,
    prob_sup_minioutcome=prob_sup_minioutcome,
    #treat_result_record=treat_result_record,
    prob_superiority=prob_superiority
    #sg= sg

  )

  return(out)

}





#' Multi.BAET: this function simulates the properties of Bayesian adaptive enrichment trial
#' @param n.sim The number of simulations
#' @param n.cores The number of cores used for computation, default is half of the total cores
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



#' @return est.mean: A ‘ns’ * ‘nt’ matrix storing the mean of the posterior mean, with row number indicating subpopulation and column number indicating treatment arm
#' @return est.sd: A ‘ns’ * ‘nt’ matrix storing the mean of the posterior SD, with row number indicating subpopulation and column number indicating treatment arm
#' @return ss.sub.trt.mean: A ‘ns’ * ‘nt’ matrix storing the averaged sample size for each subpopulation in each treatment arm, with the row number indicating the subpopulation and the column number indicating the treatment arm.
#' @return ss.sub.dist: A ‘ns’ * ‘n.sim’ matrix storing the sample size consumed for each subpopulation in each run, with row number indicating subpopulation and column number indicating the specific run
#' @return ss.sub.mean: A vector length ‘ns’ indicating the average sample size consumed in each subpopulation
#' @return ss.t.dist: A ‘n.sim’ vector storing the overall sample size consumed in each run
#' @return ss.t.mean: A number of average overall sample size consumed
#' @return power.sub: A vector length ‘ns’ indicating the (empirical) power estimated in each subpopulation
#' @return computation.time: Time consumed for computation

#' @examples
#' nt=3
#' ns=2

#' Multi.BAET(n.sim = 5,
#'            n.cores=3,
#'            nt=3, ns=2,
#'            ss.interim.es = list(c(30, 50, 100), c(20, 60, 100)),
#'            response.type = "gaussian",
#'            sig.e = 10,
#'            mean.response = matrix(c(seq(nt*ns)),nrow = nt, ncol = ns, byrow =F),
#'            prob.subpopulation = rep (1/ns, ns),
#'            prob.trtarm = rep (1/nt, nt),
#'            maxN = 100,
#'            upper = rep (1, ns),
#'            lower = rep (0, ns),
#'            rar = F,
#'            rarmin.p = 0.1,
#'            rarmax.p = 0.9,
#'            MOR = rep(0, ns),
#'            prob.MOR = rep(0.1, ns),
#'            N.MCMC = 3000,
#'            prior.cov = diag(25, ns*nt),
#'            prior.mean = rep(0, ns*nt))
#'

#' @export

Multi.BAET= function(n.sim,
                     n.cores=detectCores()/2,
                     nt=3, ns=2,
                     ss.interim.es,
                     response.type = 'gaussian',
                     sig.e = 10,
                     mean.response = matrix( c(8,12,14, 5,8,11 ), nrow = nt, ncol = ns, byrow = F),  ## input nt treatments for one subpopulation, then nt treatments for the second population
                     prob.subpopulation = c(0.6,0.4), # which population the patient belongs
                     prob.trtarm = rep (1/nt, nt),
                     maxN = 300,
                     upper = rep (0.90, ns), lower = rep (0.10, ns),
                     rar = F,
                     rarmin.p = 0.1,
                     rarmax.p = 0.9,
                     MOR = rep(-Inf, ns),
                     prob.MOR = rep(0.10, ns),
                     N.MCMC = 5000,
                     prior.cov = diag(25, ns*nt), prior.mean = c(8, 12, 14, 5, 8, 11) ){

  start = Sys.time()
  registerDoParallel(cores=n.cores)
  mt1 = foreach(i=1:n.sim, .packages = c("rjags","boot","abind","tibble", "MASS", "matlib"), .combine = "cbind") %dopar% {

    test= BAET.sim (nt=nt, ns=ns,
                    ss.interim.es = ss.interim.es,
                    response.type = response.type,
                    sig.e = sig.e,
                    mean.response = mean.response,
                    prob.subpopulation = prob.subpopulation,
                    prob.trtarm = prob.trtarm,
                    maxN = maxN,
                    upper =  upper, lower = lower,
                    rar = rar,
                    MOR = MOR,
                    prob.MOR = prob.MOR,
                    N.MCMC = N.MCMC,
                    prior.cov = prior.cov, prior.mean = prior.mean
    )


    est = matrix(rep(0, nt*ns), nrow =ns, ncol = nt)
    for (i in 1:ns){
      est[i,] = unlist(test$est[[i]])[1:nt]
    }

    sd = matrix(rep(0, nt*ns), nrow =ns, ncol = nt)
    for (i in 1:ns){
      sd[i,] = unlist(test$est[[i]])[(3*nt+1):(3*nt+nt)]
    }

    summmat = matrix(rep(0, nt*ns), nrow =ns, ncol = nt)
    summmat = test$ss.sub.trt

    ss.sub = c()
    for (i in 1:ns) {
      ss.sub[i] =sum (test$trt_sub[,2] == i)
    }
    n_terminate = test$n_terminate
    power = test$powerind

    out= list(est, sd, summmat, ss.sub, as.vector (n_terminate), as.vector(power) )
    class(out) = 'trial'
    return(out)
  }

  #treatment mean
  stacked_matrices_est <- abind(mt1[1,], along = 3)
  est = apply(stacked_matrices_est, c(1, 2), mean)
  # treament sd
  stacked_matrices_sd <- abind(mt1[2,], along = 3)
  sd = apply(stacked_matrices_sd, c(1, 2), mean)

  # numbers summary matrix
  matrix_list <- mt1[3,]  # Replace with actual matrices

  # Compute the element-wise average
  ss.sub.trt.mean <- Reduce(`+`, matrix_list) / length(matrix_list)





  # ss.sub
  ss.sub.dist = abind(mt1[4,], along = 2)
  ss.sub.mean= apply (ss.sub.dist,1,mean)
  # total.ss
  ss.t.dist = abind(mt1[5,],along = 1)
  ss.t.mean = mean(ss.t.dist)
  # power
  power.sub.dist = abind(mt1[6,], along = 2)
  power.sub= apply (power.sub.dist,1,mean)

  Sys.time()-start

  computation.time = Sys.time()-start
  out = list ( est.mean = est, est.sd =sd,
               ss.sub.trt.mean = ss.sub.trt.mean,
               ss.sub.dist = ss.sub.dist,
               ss.sub.mean=ss.sub.mean,
               ss.t.dist = ss.t.dist,
               ss.t.mean =ss.t.mean,
               power.sub = power.sub,
               computation.time = computation.time)

  return(out)

}








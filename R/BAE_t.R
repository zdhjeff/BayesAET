library(ggplot2)
library(abind)
library(reshape2)
library(rjags)
library(doParallel)
library(MASS)
library(roxygen2)
library(fastDummies)
library(parallel)


#' BAET simulator: this function simulates a Bayesian adaptive enrichment trial
#' @param nt number of treatment arms
#' @param ns number of subgroups
#' @param nb number of subjects updated at each interim look
#' @param response.type either 'binary'(probability), 'count'(lambda) or 'gaussian'
#' @param mean.response vector of mean responses: a nt (row) * ns (col) matrix, with row number indicates treatment and col number indicates subgroup
#' @param prob.subpopulation the probability of a subject coming from one specific subpopulation. default: rep (1/ns, ns)
#' @param prob.trtarm  the (initial) probability of a subject being assigned to a specific treatment arm. default: rep (1/nt, nt)
#' @param maxN the maximum sample size, trial will stop when achieving this number
#' @param upper upper probability threshold to claim superiority (in a specific subgroup); the treatment arm wins if it goes above this threshold
#' @param lower lower probability threshold to claim superiority (in a specific subgroup); the treatment arm will be dropped if it goes below this threshold
#' @param rar whether using responsive adaptive randomization (rar)
#' @param MID the minimum meaningful outcome threshold for each subgroup
#' @param prob.MID the probability threshold of being larger than the MID, treatment arms below this threshold will be dropped
#' @param N.MCMC Number of MCMC samples
#' @param prior.cov the prior covariance for a multivariate normal prior
#' @param prior.mean the prior mean for a multivariate normal prior

BAET.sim = function(nt, ns,
                    ss.interim.es = list(v1, v2),
                    response.type = "gaussian",
                    sig.e = 10,
                    mean.response = matrix(c(seq(nt*ns) ), nrow = nt, ncol = ns, byrow = F),  ## input nt treatments for one subpopulation, then nt treatments for the second population
                    prob.subpopulation = rep (1/ns, ns),
                    prob.trtarm = rep (1/nt, nt),
                    maxN = 300,
                    upper = rep (0.975, ns),
                    lower = rep (0.1, ns),
                    rar = F,
                    MID = rep(0, ns),
                    prob.MID = rep(0.5, ns),
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
  # Set a maximum sample size
  n_max <- 400
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
        treat_prob[[i]][which(treat_result[[i]]==1)]= as.vector( (apply((thetaall_new_split[[i]][,which(treat_result[[i]]==1),drop=F]-MID[i] >=0 ),2,sum))/N.MCMC ) ## minioutcome_1 is the minimum outcome needed, treatp is the probablity of larger than minioutcome_1
        treat_result[[i]]=as.numeric(treat_prob[[i]]>= prob.MID[i] ) ## whether that the treatment arm is larger than the threshold
        prob_assign[[i]][which(treat_result[[i]]==0)]=0 # this means once the trt fails in mini check, that trt stops receiving subjects

        prob_sup_minioutcome[[i]]=abind(prob_sup_minioutcome[[i]],treat_prob[[i]],along = 2) ## this is only to record treat_prob[[i]]
        treat_result_record[[i]]= abind(treat_result_record[[i]],treat_result[[i]],along = 2)## this is only to record treat_result[[i]]

        ntj[i] = length(which(apply(check[[i]], 2, mean)>lower[i] & treat_result[[i]]==1)) ## No. of groups left
        if (ntj[i]!=0){
          prob_assign[[i]] = rep(0, nt)
          if (rar == T  & ntj[[i]]>1)  {
            prob_assign[[i]][which( apply(check[[i]], 2, mean)>lower[i] & treat_result[[i]]==1)] =
              ((ntj[i] - 1)/ntj[i])*sqrt(prob_superiority[[i]][which(apply(check[[i]], 2, mean)>lower[i]& treat_result[[i]]==1),threshold_indices[i]])
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
      treat_prob[[i]][which(treat_result[[i]]==1)]= as.vector( (apply((thetaall_new_split[[i]][,which(treat_result[[i]]==1),drop=F]-MID[i] >=0 ),2,sum))/N.MCMC ) ## minioutcome_1 is the minimum outcome needed, treatp is the probablity of larger than minioutcome_1
      treat_result[[i]]=as.numeric(treat_prob[[i]]>= prob.MID[i] ) ## whether that the treatment arm is larger than the threshold
      prob_assign[[i]][which(treat_result[[i]]==0)]=0 # this means once the trt fails in mini check, that trt stops receiving subjects

      prob_sup_minioutcome[[i]]=abind(prob_sup_minioutcome[[i]],treat_prob[[i]],along = 2) ## this is only to record treat_prob[[i]]
      treat_result_record[[i]]= abind(treat_result_record[[i]],treat_result[[i]],along = 2)## this is only to record treat_result[[i]]

      ntj[i] = length(which(apply(check[[i]], 2, mean)>lower[i] & treat_result[[i]]==1)) ## No. of groups left
      if (ntj[i]!=0){
        prob_assign[[i]] = rep(0, nt)
        if (rar == T  & ntj[[i]]>1)  {
          prob_assign[[i]][which( apply(check[[i]], 2, mean)>lower[i] & treat_result[[i]]==1)] =
            ((ntj[i] - 1)/ntj[i])*sqrt(prob_superiority[[i]][which(apply(check[[i]], 2, mean)>lower[i]& treat_result[[i]]==1),threshold_indices[i]])
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

  out = list(
    n.analysis =j,
    trt_sub =trt_sub,
    est =est,
    powerind = powerind,
    y =y,
    N_terminate = length(y),
    ss.sub=counts,
    prob_sup_minioutcome=prob_sup_minioutcome,
    #treat_result_record=treat_result_record,
    prob_superiority=prob_superiority,
    sg= sg

  )

  return(out)

}

#' @return n.analysis: number of analysis conducted to end the whole trial
#' @return trt_sub: simulated treatment arm allocation (1st column) and subgroup (2nd column)
#' @return est: treatment effect estimation(posterior mean) and the 95% Credible interval bounds for each treatment arm in each subgroup
#' @return powerind: power indicator of whether the best treatment arms is correctly selected in each subgroup
#' @return y: simlulated outcome
#' @return N_terminate: the total sample size consumed when trial ends
#' @return prob_sup_minioutcome: the probability of a treatment large than the MID
#' @return prob_superiority: the probability of a treatment being the best among each subgroup



#' Multi.BAET: this function simulates the properties of Bayesian adaptive enrichment trial
#' @param n.sim number of simulations
#' @param n.cores number of cores for computation, default is half of the total cores

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
                     upper = rep (0.9, ns), lower = rep (0.0, ns),
                     rar = T,
                     MID = rep(6, ns),
                     prob.MID = rep(0.5, ns),
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
                    MID = MID,
                    prob.MID = prob.MID,
                    N.MCMC = N.MCMC,
                    prior.cov = prior.cov, prior.mean = prior.mean
    )


    nint=test$n.interim
    est = matrix(rep(0, nt*ns), nrow =ns, ncol = nt)
    for (i in 1:ns){
      est[i,] = unlist(test$est[[i]])[1:nt]
    }

    sd = matrix(rep(0, nt*ns), nrow =ns, ncol = nt)
    for (i in 1:ns){
      sd[i,] = unlist(test$est[[i]])[(3*nt+1):(3*nt+3)]
    }

    ss.sub = c()
    for (i in 1:ns) {
      ss.sub[i] =sum (test$trt_sub[,2] == i)
    }
    N_terminate = test$N_terminate
    power = test$powerind

    out= list(est, sd, ss.sub, as.vector (N_terminate), as.vector(power) )
    class(out) = 'trial'
    return(out)
  }

  #treatment mean
  stacked_matrices_est <- abind(mt1[1,], along = 3)
  est = apply(stacked_matrices_est, c(1, 2), mean)
  # treament sd
  stacked_matrices_sd <- abind(mt1[2,], along = 3)
  sd = apply(stacked_matrices_sd, c(1, 2), mean)
  # ss.sub
  ss.sub.dist = abind(mt1[3,], along = 2)
  ss.sub.mean= apply (ss.sub.dist,1,mean)
  # total.ss
  ss.t.dist = abind(mt1[4,],along = 1)
  ss.t.mean = mean(ss.t.dist)
  # power
  power.sub.dist = abind(mt1[5,], along = 2)
  power.sub= apply (power.sub.dist,1,mean)

  Sys.time()-start

  computation.time = Sys.time()-start
  out = list ( est.mean = est, est.sd =sd,
               ss.sub.dist = ss.sub.dist,
               ss.sub.mean=ss.sub.mean,
               ss.t.dist = ss.t.dist,
               ss.t.mean =ss.t.mean,
               power.sub = power.sub ,
               computation.time = computation.time)

  return(out)

}


#' @return est: the posterior mean of each treatment (col) in each subpopulation (row); e.g. est[1,2] is treatment 2 in subpopulation 1
#' @return sd: the posterior sd of each treatment (col) in each subpopulation (row)
#' @return ss.sub.dist: the sample size distribution for each subpopualtion (row)
#' @return ss.sub.mean: the expected (mean) sample size for each subpopulation
#' @return ss.t.dist: the sample size distribution for the whole trial
#' @return ss.t.mean: the mean total sample size for the whole trial
#' @return power.sub: the power for each subpopulation
#' @return computation.time: the overall simulation time used






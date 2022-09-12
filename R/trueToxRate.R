#' Simulate ground truth toxicity scenarios for BOIN dose escalation trials simulation
#'
#' Use this function to simulate toxicity-dose curves.
#' 
#' 
#' @usage trueToxRate(targetRate,ndose,sigma0,mu,sigma1,nsim,seed=NULL)
#' 
#' @param targetRate The target toxicity rate used in the BOIN design
#' @param ndose number of dose levels
#' @param sigma0 standard deviation of a normal distribution N(qnorm(\code{targetRate}),\code{sigma0}) from which the true toxicity rate of MTD is sampled. See details.
#' @param mu mean of the normal distribution N(\code{mu},\code{sigma1}^2) from which the toxicity rate difference between adjacent dose levels is sampled. See details.
#' @param sigma1 standard deviation of the normal distribution N(\code{mu},\code{sigma1}^2). See details.
#' @param nsim number of scenarios to simulate
#' @param seed the random seed for simulation

#' @details - True MTD is randomly sampled with equal probability from all dose levels, or 'all doses toxic', or 'MTD not reached'. 
#'          - The true toxicity rate of the randomly selected MTD is allowed to deviate from the target rate. This is achieved by
#'          setting p(j) = pnorm(err_j), where err_j is drawn from the normal distribution N(qnorm(\code{targetRate}),\code{sigma0}).
#'          \code{sigma0} controls the variation of the true toxicty rate of MTD from the target rate.
#'          - The difference in toxicity rate between adjacent dose levels (steepness of dose-toxicity curves) is controlled by \code{mu} and \code{sigma1}.
#'          p(j-1)=pnorm(qnorm(p(j))-err_{j-1}^2), p(j+1)=pnorm(qnorm(p(j))+err_{j+1}^2), where err_{j-1} and err_{j+1}
#'          are drawn from the normal distribution N(\code{mu},\code{sigma1}^2). 
#'          Additional constraint is applied to ensure p(j) is closest to the target rate.
#' 
#' @return \code{trueToxRate()} returns a list object, including 
#'         (1) \code{$trueRate}: a matrix of size (\code{nsim},\code{ndose}) with the true toxicity rate for each dose level for all the simulations.
#'         (2) \code{$MTD}: a vector of length \code{nsim} storing the MTD level of each simulated dose-toxicity curve. 0 if all doses are toxic, and \code{ndose+1} if MTD is not reached.
#'         (3)\code{$gap_around_target}: a vector of length \code{nsim} storing the mean toxicity rate difference around the MTD for each simulated dose-toxicity curve.
#' 
#' 
#' @author Lu, Yixing
#' 
#' @references Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I
#'             Clinical Trials, \emph{Journal of the Royal Statistical Society: Series C}, 64, 507-523.
#'
#' @export

trueToxRate<-function(targetRate,ndose,sigma0,mu,sigma1,nsim,seed=NULL){
  if (! is.null(seed)){
    set.seed(seed)
  }
  # randomly select a dose level as MTD, after introduction of a fake lowest level and a fake highest level
  fake_ndose = ndose+2
  MTD = sample(fake_ndose,nsim,replace=TRUE)
  # placeholder for true toxicity rate with fake dose levels
  trueRate  = matrix(rep(0,nsim*fake_ndose),ncol=fake_ndose)
  # placeholder for the probability difference around the target rate
  prob_gap = rep(0,nsim)
  # generate true toxicity rate for MTD for all the simulations
  err_MTD = rnorm(n=nsim,mean=qnorm(targetRate),sd=sigma0**2)
  for (i in 1:nsim){
    j = MTD[i]
    trueRate[i,j]=pnorm(err_MTD[i])
    
    if (j>1){
      for (dose in (j-1):1){
        if (dose==j-1 && err_MTD[i]>qnorm(targetRate)){
          trueRate[i,dose]=pnorm(qnorm(2*targetRate-trueRate[i,j])-rnorm(1,mean=mu,sd=sigma1**2)**2)
        }
        else{
          trueRate[i,dose]=pnorm(qnorm(trueRate[i,(dose+1)])-rnorm(1,mean=mu,sd=sigma1**2)**2)
        }
      }
    }
    if (j<fake_ndose){
      for (dose in (j+1):fake_ndose){
        if (dose==j+1 && err_MTD[i]<qnorm(targetRate)){
          trueRate[i,dose]=pnorm(qnorm(2*targetRate-trueRate[i,j])+rnorm(1,mean=mu,sd=sigma1**2)**2)
        }
        else{
          trueRate[i,dose]=pnorm(qnorm(trueRate[i,(dose-1)])+rnorm(1,mean=mu,sd=sigma1**2)**2)
        }
      }
    }
    # calculate the average probability difference around the target rate: (|p(j-1)-p(j)|+|p(j+1)-p(j)|)/2
    diff_lower = ifelse(j>1,abs(trueRate[i,j-1]-trueRate[i,j]),0)
    diff_upper = ifelse(j<fake_ndose,abs(trueRate[i,j+1]-trueRate[i,j]),0)
    prob_gap[i]= mean(diff_lower,diff_upper)
  }
  
  # output
  trueRate = as.data.frame(trueRate[,2:(fake_ndose-1)])
  colnames(trueRate)<-paste('Dose',1:ncol(trueRate))
  MTD=MTD-1
  #MTD[which(MTD==0)]='All doses toxic'
  #MTD[which(MTD==ndose+1)]='MTD not reached'
  avg_prob_gap = mean(prob_gap)
  out=list(trueRate=trueRate,MTD=MTD,gap_around_target=avg_prob_gap)
  return(out)
}

#' 
#' Conduct gridsearch for BOIN dose escalation trials with specified parameters
#' 
#' Apply \code{run_boinsim} function to all combinations of parameters specified in the argument. 
#'
#' @usage gridsearch_boinsim(object,trueRate,
#'                           cohortSize,min_cohortEnrolled,min_cohortEval,maxSubjects,
#'                           dropout_rate,sf_rate,accrual_rate,
#'                           window_screen,window_DLT,meanDLTtime,lagPlanned,lagDropout,
#'                           k.earlystop,n.earlystop,ncores=1)
#' @param object a boinsim object that is output by the function run_boinsim
#' @param trueRate vector of true DLT rates for each dose level. 
#' @param cohortSize number of patients recruited at a time during the dose escalation trial
#' @param min_cohortEnrolled minimum number of patients passing the screening at each cohort
#' @param min_cohortEval minimum number of patients needed at a cohort if dropouts occur
#' @param maxSubjects If this many subjects completed the DLT assessment, the end of escalation is triggered.
#' @param accrual_rate the rate at which new patients are recruited (per day)
#' @param sf_rate probability of screen failure
#' @param dropout_rate probability of dropout 
#' @param window_screen the screening window (days)
#' @param window_DLT the assessment window for DLTs (days)
#' @param meanDLTtime mean time for DLT development (among patients who do develop DLT reaction)
#' @param lagPlanned the delay or advance (in days) of recruitment for the next cohort from the completion of assessment 
#' @param lagDropout the delay (in days) of searching for replacement after a dropout is confirmed
#' @param k.earlystop can be FALSE or a non-negative integer. When it is a non-negative integer, it controls the maximum number of cohorts enrolled at each dose level. One of k.earlystop and n.earlystop needs to be FALSE. 
#' @param n.earlystop can be a non-negative integer or FALSE. When not FALSE, it controls the number of subjects enrolled in one dose level which triggers end of escalation.
#' @param ncores Number of cores available to run the simulations in parallel.Default to 1, when the jobs will be run in sequence. 
#' 
#' @details The parameters that take a single value in \code{run_boinsim} could be provided as a single value or a vector containing all the values to be screened.
#'          If more than one set of \code{trueRate} is needed, it could be provided as a data frame with column names: Dose 1, ... Dose k,(index); or a (named) list of vectors.
#'          The "index" column or the names of the list is optional. If available, it will be used as an unique identifier of the true toxicity scenario. 
#'          
#'          All unique combinations of the following parameters will be included in the searching grid: trueRate,
#'          nesting(cohortSize,min_cohortEnrolled,min_cohortEval),maxSubjects,
#'          dropout_rate,sf_rate,accrual_rate,window_screen,window_DLT,meanDLTtime,lagPlanned,lagDropout,nesting(k.earlystop,n.earlystop). 
#'          If any of the above parameters is not specified in the argument, the parameter value from the example boinsim \code{object} will be used.
#'        
#'          - Note: for the following two groups of parameters, only combinations already present in the argument (parameter values from the same vector index position) will be included:
#'          (1) \code{cohortSize},\code{min_cohortEnrolled},and\code{min_cohortEval}
#'          (2) \code{k.earlystop} and \code{n.earlystop}
#'          
#' @return a tibble object with each row storing the boinsim output for a unique parameter combination
#' @import dplyr tidyr future furrr parallel foreach doParallel
#' @export

gridsearch_boinsim<-function(object,
                             trueRate,
                             cohortSize,min_cohortEnrolled,min_cohortEval,maxSubjects,
                             dropout_rate,sf_rate,accrual_rate,
                             window_screen,window_DLT,meanDLTtime,lagPlanned,lagDropout,
                             k.earlystop,n.earlystop,ncores=1){
  # simulation setup 
  sim.setup = object$sim.setup[[1]]
  # the following parameters can be included in the searching grid
  if (missing(trueRate)) trueRate = list(sim.setup$trueRate)
  # when trueRate is input as a data frame
  if (inherits(trueRate,'data.frame')){
    if (is.null(trueRate$index)){index = rownames(trueRate)}
    else {index = trueRate$index}
    trueRate = trueRate%>%select(starts_with('Dose'))%>%t()%>%as.data.frame()%>%as.list()
  }
  # when trueRate is input as a list
  else if (inherits(trueRate,'list')){
    if (is.null(names(trueRate))){index=1:length(trueRate)}
    else{index=names(trueRate)}
  }
  else{stop('trueRate needs to be either a list or a data.frame.')}
  
  
  if (missing(cohortSize)) cohortSize = unique(sim.setup$cohortSize)
  if (missing(min_cohortEnrolled)) min_cohortEnrolled=unique(sim.setup$min_cohortEnrolled)
  if (missing(min_cohortEval)) min_cohortEval=unique(sim.setup$min_cohortEval)
  if (missing(maxSubjects)) maxSubjects=unique(sim.setup$maxSubjects)
  if (missing(dropout_rate)) dropout_rate=unique(sim.setup$dropout_rate)
  if (missing(sf_rate)) sf_rate=unique(sim.setup$sf_rate)
  if (missing(accrual_rate)) accrual_rate=unique(sim.setup$accrual_rate)
  if (missing(window_screen)) window_screen=unique(sim.setup$window_screen)
  if (missing(window_DLT)) window_DLT=unique(sim.setup$window_DLT)
  if (missing(meanDLTtime)) meanDLTtime=unique(sim.setup$meanDLTtime)
  if (missing(lagPlanned)) lagPlanned=unique(sim.setup$lagPlanned)
  if (missing(lagDropout)) lagDropout=unique(sim.setup$lagDropout)
  if (missing(k.earlystop)) k.earlystop=unique(sim.setup$k.earlystop)
  if (missing(n.earlystop)) n.earlystop=unique(sim.setup$n.earlystop)
  # the following parameters will be fixed and inherited from the example boinsim object
  ntrial=unique(sim.setup$ntrial)
  targetRate = unique(sim.setup$targetRate)
  p.saf=unique(sim.setup$p.saf)
  p.tox=unique(sim.setup$p.tox)
  cutoff.eli=unique(sim.setup$cutoff.eli)
  seed=unique(sim.setup$seed)
  # construct a parameter grid 
  parameter_grid = crossing(nesting(trueRate,index),
                            nesting(cohortSize,min_cohortEnrolled,min_cohortEval),maxSubjects,
                            dropout_rate,sf_rate,accrual_rate,window_screen,window_DLT,meanDLTtime,lagPlanned,lagDropout,
                            nesting(k.earlystop,n.earlystop))
  # map run_boinsim to each row of parameter_grid
  if (ncores==1){
    Raw<-map(.x=1:dim(parameter_grid)[1],
             ~c(run_boinsim(targetRate=targetRate,trueRate=parameter_grid$trueRate[[.x]], 
                            cohortSize=parameter_grid$cohortSize[[.x]],min_cohortEnrolled=parameter_grid$min_cohortEnrolled[[.x]],
                            min_cohortEval=parameter_grid$min_cohortEval[[.x]],maxSubjects=parameter_grid$maxSubjects[[.x]],
                            dropout_rate=parameter_grid$dropout_rate[[.x]],sf_rate=parameter_grid$sf_rate[[.x]],accrual_rate=parameter_grid$accrual_rate[[.x]],
                            window_screen=parameter_grid$window_screen[[.x]],window_DLT=parameter_grid$window_DLT[[.x]],meanDLTtime=parameter_grid$meanDLTtime[[.x]],
                            lagPlanned=parameter_grid$lagPlanned[[.x]],lagDropout=parameter_grid$lagDropout[[.x]],
                            k.earlystop=parameter_grid$k.earlystop[[.x]],n.earlystop=parameter_grid$n.earlystop[[.x]],
                            p.saf=p.saf,p.tox=p.tox,cutoff.eli=cutoff.eli,ntrial=ntrial,seed=seed,rawOutput=FALSE),
                scenario_index=parameter_grid$index[.x]),.progress=TRUE)
    Raw<-Raw%>%bind_rows()
  }
  else {
    cl<-parallel::makeCluster(ncores)
    future::plan(strategy='cluster',workers=cl)
    options(future.rng.onMisuse='ignore')
    Raw<-furrr::future_map(.x=1:dim(parameter_grid)[1],
             ~c(run_boinsim(targetRate=targetRate,trueRate=parameter_grid$trueRate[[.x]], 
                            cohortSize=parameter_grid$cohortSize[.x],min_cohortEnrolled=parameter_grid$min_cohortEnrolled[.x],
                            min_cohortEval=parameter_grid$min_cohortEval[.x],maxSubjects=parameter_grid$maxSubjects[.x],
                            dropout_rate=parameter_grid$dropout_rate[.x],sf_rate=parameter_grid$sf_rate[.x],accrual_rate=parameter_grid$accrual_rate[.x],
                            window_screen=parameter_grid$window_screen[.x],window_DLT=parameter_grid$window_DLT[.x],meanDLTtime=parameter_grid$meanDLTtime[.x],
                            lagPlanned=parameter_grid$lagPlanned[.x],lagDropout=parameter_grid$lagDropout[.x],
                            k.earlystop=parameter_grid$k.earlystop[.x],n.earlystop=parameter_grid$n.earlystop[.x],
                            p.saf=p.saf,p.tox=p.tox,cutoff.eli=cutoff.eli,ntrial=ntrial,seed=seed,rawOutput=FALSE),
                scenario_index=parameter_grid$index[.x]),.progress=TRUE)
    splits = parallel::splitIndices(length(Raw),ncores) 
    doParallel::registerDoParallel(cl)
    Raw = foreach::foreach(i = 1:length(splits)) %dopar% {
      bind_rows(Raw[splits[[i]]])
    }
    Raw=bind_rows(Raw)
    parallel::stopCluster(cl)
  }
  return(Raw)
}
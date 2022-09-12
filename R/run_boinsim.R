#'
#' Simulate BOIN dose escalation with patient dropout and timeline monitoring
#'
#' @usage run_boinsim(targetRate, trueRate, cohortSize, min_cohortEnrolled,
#'        min_cohortEval,maxSubjects,dropout_rate,sf_rate=0.2,accrual_rate=0.1,
#'        window_screen=14,window_DLT=28,meanDLTtime=7,
#'        lagPlanned=0,lagDropout=0,k.earlystop=3,n.earlystop=FALSE,
#'        p.saf=0.6*targetRate,p.tox=1.227*targetRate,cutoff.eli=0.95,
#'        ntrial=100,seed=2022,rawOutput=FALSE)

#' @param targetRate The target rate for the BOIN design
#' @param trueRate vector of true DLT rates for each dose level
#' @param cohortSize number of patients recruited at a time during the dose escalation trial
#' @param min_cohortEnrolled minimum number of patients passing the screening at each cohort
#' @param min_cohortEval minimum number of patients needed at a cohort if dropouts occur
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
#' @param maxSubjects If this many subjects completed the DLT assessment, the end of escalation is triggered.
#' @param p.saf the highest toxicity rate that is deemed subtherapeutic (lower boundary of the target interval)
#' @param p.tox the lowest toxicity rate that is deemed overly toxic (upper boundary of the target interval)
#' @param cutoff.eli The elimination threshold used to generate the BOIN decision table. Default to 0.95
#' @param ntrial number of simulation trials to run. Default to 100.
#' @param seed the random seed for simulation. Defaul to 2022.
#' @param rawOutput whether to return the raw output. Default to FALSE, when only a summary of all simulation trials is returned.

#' @return If rawOutput==TRUE, the output for a single simulation trial is returned in a list, including:
#'         (1) DLTS: number of patients experiencing DLT at each dose level
#'         (2) DLTs_interval: number of patients experiencing DLT at each interval level (below, within, or above the target interval)
#'         (3) Subjects/Subjects_perc: number/percentage of patients recruited at each dose level
#'         (4) SubjectsEnrolled/SubjectsEnrolled_perc: number/percentage of patients who pass the screening at each dose level
#'         (5) SubjectsCompleted/SubjectsCompleted_perc: number/percentage of evaluable patients at each dose level
#'         (6) SubjectsInterval/SubjectsInterval_perc: number/percentage of patients recruited at each interval level
#'         (7) SubjectsEnrolledInterval/SubjectsEnrolledInterval_perc: number/percentage of patients who pass the screening at each interval level
#'         (8) SubjectsCompletedInterval/SubjectsCompletedInterval_perc: number/percentage of evaluable patients at each interval level
#'         (9) Cohorts: number of cohorts assigned at each dose level
#'         (10) CohortsInterval: number of cohorts assigned at each interval level
#'         (11) avg_DurationCohort: average duration for a cohort to complete, from searching patients to completion of assessment
#'         (12) DurationTrial: total duration of the trial
#'         (13) DosePath: the path of dose escalation/deescaltion/stay during the trial
#'         (14) DoseSelected: dose selected as MTD
#'         (15) DoseSelected_interval: whether the selected MTD is below, within, or above the target interval
#'         
#'         If rawOutput==FALSE (default), a mean over all simulation trials are calculated and output for item 1~12, in addition the following items are included in the output list:
#'         (13) SlcPerc_dose: MTD selection percentage for 'All excluded' and all dose levels
#'         (14) SlcPerc_interval: MTD selection percentage for 'All excluded' and below, within, or above the target interval
#'         
#'         A dataframe of simulation setup (sim.setup) will be included in the output list for both versions.
#'         
#' 
#' @author Lu, Yixing; Henner, William
#' 
#' @import BOIN
#' @import dplyr
#' 
#' @export

run_boinsim<-function(targetRate, trueRate, cohortSize,min_cohortEnrolled,min_cohortEval,maxSubjects,dropout_rate,sf_rate=0.2,
         accrual_rate=0.1,window_screen=14,window_DLT=28,meanDLTtime=7,lagPlanned=0,lagDropout=0,
         k.earlystop=3,n.earlystop=FALSE,p.saf=0.6*targetRate,p.tox=1.227*targetRate,cutoff.eli=0.95,ntrial=100,seed=2022,rawOutput=FALSE){
  
  # sanity check 
  if (cohortSize<min_cohortEnrolled|cohortSize<min_cohortEval|min_cohortEnrolled<min_cohortEval){
    stop("The inequality: cohortSize>=min_cohortEnrolled>=min_cohortEval needs to be fulfilled.")
  }
  if (! is.null(seed)){
    set.seed(seed)
  }
  # Derive BOIN decision table
  DecisionTable=BOIN::get.boundary(targetRate,100,cohortsize=3,cutoff.eli = cutoff.eli,p.saf=p.saf,p.tox=p.tox)
  escB = DecisionTable$lambda_e
  deB = DecisionTable$lambda_d
  ndose = length(trueRate)
  
  # set up list to store simulation output for ntrial runs
  rawOut = list()
  
  # simulation setup
  sim.setup = tibble(targetRate=targetRate,
                     trueRate=list(as.numeric(trueRate)),
                     cohortSize=cohortSize,
                     min_cohortEnrolled=min_cohortEnrolled,
                     min_cohortEval=min_cohortEval,
                     maxSubjects=maxSubjects,
                     dropout_rate=dropout_rate,
                     sf_rate=sf_rate,
                     accrual_rate=accrual_rate,
                     window_screen=window_screen,
                     window_DLT=window_DLT,
                     meanDLTtime=meanDLTtime,
                     lagPlanned=lagPlanned,
                     lagDropout=lagDropout,
                     k.earlystop=k.earlystop,
                     n.earlystop=FALSE,
                     cutoff.eli=cutoff.eli,
                     p.saf=p.saf,
                     p.tox=p.tox,
                     ntrial=ntrial,
                     seed=seed)
  
  ## helper function: summarize dose-wise outputs into interval-wise ones
  assign_interval <- function(output,TargetInterval){
    outputInterval = rep(0,3)
    outputInterval[1]=sum(output[TargetInterval==0]) # output at doses below target interval
    outputInterval[2]=sum(output[TargetInterval==1]) # output at doses within target interval
    outputInterval[3]=sum(output[TargetInterval==2]) # output at doses above target interval
    names(outputInterval)=c('Below target interval','Within target interval','Above target interval')
    return(outputInterval)
  }
  
  for (trial in 1:ntrial){
    # randomly generate screen failure event (=1 if failed)
    random_sf = rbinom(maxSubjects,1,sf_rate)
    # randomly generate dropout events (=1 if passes the screening but drops out, =0 otherwise)
    random_dropouts=rep(0,maxSubjects)
    random_dropouts[random_sf==0]=rbinom(sum(random_sf==0),1,dropout_rate)
    while (sum(1-(random_sf+random_dropouts))<maxSubjects) {
      new=rbinom(1,1,sf_rate)
      random_sf = c(random_sf,new)
      random_dropouts = c(random_dropouts,ifelse(new==1,0,rbinom(1,1,dropout_rate)))
    } # total evaluable patients = maxSubjects
    totSubjects = length(random_dropouts)
    
    # randomly generate search time for each patient
    random_searchTime = ceiling(rexp(totSubjects,accrual_rate))  
    # number of days during screening until failure(=window_screen if pass the screening)
    random_sfDay=sample(window_screen,totSubjects,replace=TRUE)
    
    # number of days to dropout (=0 if do not drop out)
    random_dropoutDay = rep(0,totSubjects)
    random_dropoutDay[which(random_dropouts==1)]=sample(window_DLT-1,sum(random_dropouts==1),replace = TRUE)
    
    # randomly generate DLT results
    random_DLT = runif(totSubjects)
    random_DLT[which(random_sf==1 | random_dropouts==1)]=1 
    
    # assert one and only one of k.earlystop and n.earlystop is FALSE
    stopifnot((k.earlystop!=FALSE && n.earlystop==FALSE)||(k.earlystop==FALSE && n.earlystop!=FALSE))
    
    # Set up tracking vectors
    I = 1 # subject number tracker
    DLTs = rep(0,ndose) # Total number of DLTs at each dose level
    dose = 1 # Current dose level
    DoseElim = rep(F,ndose) # Eliminated dose levels
    Subjects = rep(0,ndose) # Total number of subjects recruited at each dose level
    SubjectsEnrolled = rep(0,ndose) # Total number of subjects passing screening and enrolled for DLT assessment
    SubjectsCompleted = rep(0,ndose) # Total number of subjects completed the DLT assessment at each dose level
    Cohorts = rep(0,ndose) # Total number of cohorts enrolled at each dose level
    Status_DLT=rep(0,totSubjects) # Whether each patient develop DLT or not
    
    size = cohortSize # number of subjects to recruit in a cohort
    
    DLTrate = NA # DLT rate in current cohort
    
    DurationCohort  = NULL # duration of each cohort
    
    DurationSubject = random_searchTime+random_sfDay+random_dropoutDay # duration of each patient
    
    #StatusSubject = random_dropouts # status code of patients: 0=cleared DLT, 1=dropouts, 2=DLT
    
    timeStamp2 = NULL # time stamp used to calculate total trial duration
    
    DosePath = c(dose) # track the dose level path during dose escalation
    
    while(TRUE){
      #The numbers of the subjects in this cohort
      if (n.earlystop!=FALSE){
        cohortSizeTrue = min(size,totSubjects-I+1,n.earlystop-SubjectsCompleted[dose])
      }
      else if (k.earlystop!=FALSE){
        cohortSizeTrue = min(size,totSubjects-I+1)
      }
      # checkpoint: if number of available patients at recruitment is smaller than the min. evaluable patients required per cohort, break
      if (cohortSizeTrue<min_cohortEval){
        break
      }
      else{
        sI <- I:(I+cohortSizeTrue-1)
      }
      
      # search time+screen failure day (window_screen if passes screening) for patients in the current cohort
      timeStamp1_cohort=random_searchTime[sI]+random_sfDay[sI]
      # update the trackers of number of patients at the current dose level
      Subjects[dose]=Subjects[dose]+cohortSizeTrue
      Cohorts[dose]=Cohorts[dose]+1
      # update the patient number tracker
      I=max(sI)+1
      
      # check if number of patients passed the screening in the current cohort drops below min_cohortEnrolled and needs replacement
      while (sum(random_sf[sI]==0)<min_cohortEnrolled
             & any(random_sf[sI]==1)
             & I <= totSubjects) {
        # determine the patient ID of the one after whose failure we will need to search for replacement
        failures = which(random_sf[sI]==1)
        i=failures[sort(timeStamp1_cohort[failures],decreasing=TRUE,index.return=TRUE)$ix[min_cohortEnrolled-(cohortSizeTrue-length(failures))]]
        sI[i]=I 
        # update the cohort duration for this slot, including a lag time, and search time for a replacement
        timeStamp1_cohort[i]=timeStamp1_cohort[i]+lagDropout+random_searchTime[I]+random_sfDay[I]
        # update total patients recruited for the current dose
        Subjects[dose]=Subjects[dose]+1
        # update subject number tracker
        I=I+1
      }
      
      # select for patients who passed the screening in this cohort
      enrolled=(random_sf[sI]==0)
      # update tracker of number and ID of enrolled patients at the current dose level who pass the screening
      SubjectsEnrolled[dose]=SubjectsEnrolled[dose]+sum(enrolled)
      sI = sI[enrolled]
      
      # only select those who pass the screening (same length as current sI)
      timeStamp1_enrolled=timeStamp1_cohort[enrolled]+random_dropoutDay[sI]
      
      # check if number of enrolled patients in the current cohort drops below min_cohortEval and needs replacement
      while (sum(random_dropouts[sI]==0)<min_cohortEval
             & any(random_dropouts[sI]==1)
             & I <= totSubjects) {
        dropouts = which(random_dropouts[sI]==1)
        i=dropouts[sort(timeStamp1_enrolled[dropouts],decreasing=TRUE,index.return=TRUE)$ix[min_cohortEval-(sum(enrolled)-length(dropouts))]]
        
        # search until finding a replacement who passes the screening; update the cohort duration for this slot, including a lag time, search time for a replacement, and screening
        while (random_sf[I]==1){
          timeStamp1_enrolled[i]=timeStamp1_enrolled[i]+lagDropout+random_searchTime[I]+random_sfDay[I]
          # update total recruited patients for the current dose
          Subjects[dose]=Subjects[dose]+1
          I=I+1
        }
        # refill the ith slot in this cohort with a patient passing the screening; update cohort duration with dropout day
        timeStamp1_enrolled[i]=timeStamp1_enrolled[i]+lagDropout+random_searchTime[I]+random_sfDay[I]+random_dropoutDay[I]
        sI[i]=I # enroll a new patient who passed the screening as replacement
        # update the tracker for number of recruited and enrolled subjects at the current cohort
        Subjects[dose]=Subjects[dose]+1
        SubjectsEnrolled[dose]=SubjectsEnrolled[dose]+1
        # update subject ID tracker
        I=I+1
      }
      # update timeStamp1_cohort to include the dropoutDay for those passed the screening
      timeStamp1_cohort[enrolled]=timeStamp1_enrolled
      
      # select evaluable patients in this cohort
      evaluable = (random_dropouts[sI]==0)
      # update tracker of number and ID of evaluable patients at the current dose level 
      SubjectsCompleted[dose] = SubjectsCompleted[dose]+sum(evaluable)
      sI=sI[evaluable]
      
      # DLT results for all patients completed the assessment in this cohort (dropouts would be counted as non-DLT)
      ss = random_DLT[sI]<as.numeric(trueRate[dose])
      dayToDLT = rep(window_DLT,length(ss))
      tmp = ceiling(rexp(1,1/meanDLTtime))
      dayToDLT[ss]=ifelse(tmp<=window_DLT,tmp,window_DLT)
      
      # update subject DLT status
      Status_DLT[sI][ss==1]=1
      expectedEnd = timeStamp1_cohort+window_DLT+lagPlanned # expected end date of the cohort: time to enroll all evaluable patients+window_DLT+lagPlanned
      # update cohort duration 
      timeStamp1_cohort[enrolled][evaluable]=timeStamp1_cohort[enrolled][evaluable]+dayToDLT
      DurationCohort = c(DurationCohort,max(timeStamp1_cohort))
      # update subject duration
      DurationSubject[sI]=DurationSubject[sI]+dayToDLT
      # time point to start searching for the next cohort
      trueEnd = timeStamp1_cohort+max(0,lagPlanned) # true end date the of cohort
      timeStamp2 = c(timeStamp2,max(apply(cbind(expectedEnd,trueEnd),1,min))) # start searching for the next cohort whichever comes first
      
      #add the DLTs to the tracking vector
      DLTs[dose] = DLTs[dose] + sum(ss)
      
      # DLT rate in current cohort: terminate if total evaluable patients at the dose level<min_cohortEval
      if (SubjectsCompleted[dose]>=min_cohortEval){
        DLTrate = DLTs[dose]/SubjectsCompleted[dose]
      }
      else{
        break
      }
      
      #Check Deescalate
      if(DLTrate >= deB){
        #check Elimination (Assumes elimination limit is higher than Deescalate)
        if(1-pbeta(q = targetRate,1+DLTs[dose],1+SubjectsCompleted[dose]-DLTs[dose]) >cutoff.eli & SubjectsCompleted[dose] >2){
          DoseElim[dose:ndose] = T
          # End criteria: Dose elimination when at the lowest dose level.
          if(dose== 1){
            #print('Escalation terminated: lowest dose eliminated')
            break
          }
        }
        else{
          dose = max(dose -1,1) # if already at the lowest dose, will enroll more at lowest dose level instead of deescalating
        }
      } 
      
      #Check for escalation
      else if (DLTrate <= escB) {
        if (dose>=max(which(!DoseElim)) 
        ){
          #print('Escalation terminated: earlystop criteria is triggered')
          break
        }
        else{
          dose = min(dose+1,max(which(!DoseElim))) # will enroll more at the un-eliminated highest level instead of escalating
        }
      }
      
      # check for stay
      else{
        dose=dose
      }
      # update dose tracker
      DosePath = c(DosePath,dose)
      
      # End criteria: Max subjects exceeded overall
      if(I > totSubjects | (SubjectsCompleted[dose]>=n.earlystop & Cohorts[dose]>=k.earlystop)){
        #print('Escalation terminated: exceeding maximum subjects')
        break
      }
    }
    
    # Select the highest dose level which is below the upper interval bound as MTD.
    # Visited = which(SubjectsCompleted >=6) # AND at least 6 patients have completed the DLT assessment at the selected level
    # Selection <- Visited[which(DLTs[Visited]/SubjectsCompleted[Visited]<deB)] 
    # Selection <- ifelse(length(Selection) >0,max(Selection),0)
    
    # select MTD using isotonic regression with the constraint that the isotonic estimate of toxicity rate must be lower than the deescalation boundary
    Selection <- BOIN::select.mtd(targetRate,SubjectsCompleted,DLTs,boundMTD = TRUE,p.tox=p.tox)
    Selection<-Selection$MTD
    Selection=ifelse(Selection==99,0,Selection)
    #---------------------------------------------------------------------------------------------------------------------
    # prepare outputs
    # assign each dose level to below, within, or above the target interval
    TargetInterval=(as.numeric(trueRate)>p.saf)+(as.numeric(trueRate)>p.tox)
    
    # assign selected dose level to below, within, or above the target interval
    Selection_interval=ifelse(Selection==0, 0, TargetInterval[Selection]+1)
    Selection_interval=factor(Selection_interval,levels=c(0,1,2,3),labels=c('All excluded','Below target interval','Within target interval','Above target interval'))
    
    # calculate total trial duration
    DurationTrial = sum(timeStamp2[-length(timeStamp2)])+DurationCohort[length(DurationCohort)]
    
    # total DLTs, Subjects, SubjectsEnrolled,SubjectsCompleted, cohorts
    names(DLTs)=paste('Dose',1:ndose)
    names(Subjects)=paste('Dose',1:ndose)
    names(SubjectsEnrolled)=paste('Dose',1:ndose)
    names(SubjectsCompleted)=paste('Dose',1:ndose)
    names(Cohorts)=paste('Dose',1:ndose)
    # percentage of subjects assigned to each dose
    Subjects_perc=Subjects/sum(Subjects)*100
    SubjectsEnrolled_perc=SubjectsEnrolled/sum(SubjectsEnrolled)*100
    SubjectsCompleted_perc=SubjectsCompleted/sum(SubjectsCompleted)*100
    
    # calculate interval-based Subjects, SubjectsEnrolled,SubjectsCompleted,Cohorts, DLTs
    DLTs_interval=assign_interval(DLTs,TargetInterval)
    SubjectsInterval=assign_interval(Subjects,TargetInterval)
    SubjectsEnrolledInterval=assign_interval(SubjectsEnrolled,TargetInterval)
    SubjectsCompletedInterval=assign_interval(SubjectsCompleted,TargetInterval)
    SubjectsInterval_perc=assign_interval(Subjects_perc,TargetInterval)
    SubjectsEnrolledInterval_perc=assign_interval(SubjectsEnrolled_perc,TargetInterval)
    SubjectsCompletedInterval_perc=assign_interval(SubjectsCompleted_perc,TargetInterval)
    CohortsInterval=assign_interval(Cohorts,TargetInterval)
    
    # final output for one simulation trial
    out=list(sim.setup=list(sim.setup),DLTs = list(DLTs), DLTs_interval=list(DLTs_interval),
             Subjects=list(Subjects), Subjects_perc=list(Subjects_perc),SubjectsInterval=list(SubjectsInterval),SubjectsInterval_perc=list(SubjectsInterval_perc),
             SubjectsCompleted=list(SubjectsCompleted),SubjectsCompleted_perc=list(SubjectsCompleted_perc),SubjectsCompletedInterval=list(SubjectsCompletedInterval),SubjectsCompletedInterval_perc=list(SubjectsCompletedInterval_perc),
             SubjectsEnrolled=list(SubjectsEnrolled),SubjectsEnrolled_perc=list(SubjectsEnrolled_perc),SubjectsEnrolledInterval=list(SubjectsEnrolledInterval),SubjectsEnrolledInterval_perc=list(SubjectsEnrolledInterval_perc),
             Cohorts=list(Cohorts), CohortsInterval=list(CohortsInterval),
             DosePath = list(DosePath),
             avg_DurationCohort=mean(DurationCohort),DurationTrial=DurationTrial, DoseSelected = Selection, DoseSelected_interval=Selection_interval)
    # append to the list of all trials' output
    rawOut[[trial]]=out
  }
  rawOut <- rawOut%>%bind_rows()
  
  # summarize output
  # % MTD selection
  SlcPerc_dose = rep(0,ndose+1)
  for (i in 0:ndose){
    SlcPerc_dose[i+1]=sum(rawOut$DoseSelected==i)/ntrial*100
  }
  names(SlcPerc_dose)=c('All excluded',paste('Dose',1:ndose))
  
  intervalName <- c('All excluded','Below target interval','Within target interval','Above target interval') 
  SlcPerc_interval = rep(0,length(intervalName))
  for (i in 1:length(intervalName)){
    SlcPerc_interval[i]=sum(rawOut$DoseSelected_interval==intervalName[i])/ntrial*100
  }
  names(SlcPerc_interval)=intervalName
  
  ## helper function: summarize output dose-wise or interval-wise
  levelwise.summary<-function(Raw,colname,type,ndose){
    DoseName = paste('Dose',1:ndose)
    intervalName <- c('Below target interval','Within target interval','Above target interval')
    n = do.call(rbind,Raw%>%pull(colname))
    n = colMeans(cbind(n,apply(n,1,sum))) 
    if (type=='dose'){
      names(n)=c(DoseName,'Total')
    }
    else if (type=='interval'){
      names(n)=c(intervalName,'Total')
    }
    return(list(n))
  }
  
  summaryOut=list(sim.setup=list(sim.setup), 
           DLTs = levelwise.summary(rawOut,'DLTs',type='dose',ndose=ndose), DLTs_interval=levelwise.summary(rawOut,'DLTs_interval',type='interval',ndose=ndose),
           Subjects=levelwise.summary(rawOut,'Subjects',type='dose',ndose=ndose), Subjects_perc=levelwise.summary(rawOut,'Subjects_perc',type='dose',ndose=ndose),SubjectsInterval=levelwise.summary(rawOut,'SubjectsInterval',type='interval',ndose=ndose),SubjectsInterval_perc=levelwise.summary(rawOut,'SubjectsInterval_perc',type='interval',ndose=ndose),
           SubjectsCompleted=levelwise.summary(rawOut,'SubjectsCompleted',type='dose',ndose=ndose),SubjectsCompleted_perc=levelwise.summary(rawOut,'SubjectsCompleted_perc',type='dose',ndose=ndose),SubjectsCompletedInterval=levelwise.summary(rawOut,'SubjectsCompletedInterval',type='interval',ndose=ndose),SubjectsCompletedInterval_perc=levelwise.summary(rawOut,'SubjectsCompletedInterval_perc',type='interval',ndose=ndose),
           SubjectsEnrolled=levelwise.summary(rawOut,'SubjectsEnrolled',type='dose',ndose=ndose),SubjectsEnrolled_perc=levelwise.summary(rawOut,'SubjectsEnrolled_perc',type='dose',ndose=ndose),SubjectsEnrolledInterval=levelwise.summary(rawOut,'SubjectsEnrolledInterval',type='interval',ndose=ndose),SubjectsEnrolledInterval_perc=levelwise.summary(rawOut,'SubjectsEnrolledInterval_perc',type='interval',ndose=ndose),
           Cohorts=levelwise.summary(rawOut,'Cohorts',type='dose',ndose=ndose), CohortsInterval=levelwise.summary(rawOut,'CohortsInterval',type='interval',ndose=ndose),
           SlcPerc_dose = list(SlcPerc_dose), SlcPerc_interval=list(SlcPerc_interval),
           avg_DurationCohort=mean(rawOut$avg_DurationCohort),DurationTrial=mean(rawOut$DurationTrial))
  if (rawOutput){out=rawOut}
  else {out=as_tibble(summaryOut)}
  return(out)
}
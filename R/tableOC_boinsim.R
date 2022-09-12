#' Generate operating characteristics summary tables.
#' 
#' Generate a summary table with true toxicity rates, number of DLTs, number of recruited/enrolled/evaluable subjects, and MTD selection percentage at each dose level.
#' 
#' @usage tableOC_boinsim(object,
#'               cohortSize,min_cohortEnrolled,min_cohortEval,maxSubjects,
#'               dropout_rate,sf_rate,accrual_rate,
#'               window_screen,window_DLT,meanDLTtime,lagPlanned,lagDropout,
#'               k.earlystop,n.earlystop,tableOutput=TRUE)
#' @param object an object output by run_boinsim or gridsearch_boinsim. 
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
#' @param tableOutput boolean. Default to TRUE, when formatted tables will be output. If FALSE, a list of data frames will be output.
#' @details - If \code{object} is an output by run_boinsim, \code{rawOutput} has to be set to TRUE in run_boinsim. The plots will be a summary of all simulation trials.
#'            There is no need to specify parameters other than \code{type} and \code{which}.
#'          - If \code{object} is an output by gridsearch_boinsim, the plots will be a summary across all toxicity scenarios, but only for one combination of other parameters as specified by the rest of the arguments.
#'            For the arguments that only take one value in \code{object} (not subject to gridsearch), there is no need to specify here again. 
#' 
#' @import dplyr tidyr flextable officer caret
#' @export

tableOC_boinsim<-function(object,
                           cohortSize,min_cohortEnrolled,min_cohortEval,maxSubjects,
                           dropout_rate,sf_rate,accrual_rate,
                           window_screen,window_DLT,meanDLTtime,lagPlanned,lagDropout,
                           k.earlystop,n.earlystop,tableOutput=TRUE){
  
  # table format settings
  flextable::set_flextable_defaults(
    font.family = "Times New Roman", 
    font.size = 10,
    font.color = "black",
    table.layout = "fixed"
  )
  
  # filter dataset with specified parameters
  sim.setup=object$sim.setup%>%bind_rows()
  ## helper function
  unique_check<-function(sim.setup,parameter){
    tmp=sim.setup%>%pull(parameter)
    tmp=unique(tmp)
    if (length(tmp)!=1){stop('Only one value is allowed for each parameter.')}
    return(tmp)
  }
  if (missing(cohortSize)) cohortSize=unique_check(sim.setup,'cohortSize')
  if (missing(min_cohortEnrolled)) min_cohortEnrolled=unique_check(sim.setup,'min_cohortEnrolled')
  if (missing(min_cohortEval)) min_cohortEval=unique_check(sim.setup,'min_cohortEval')
  if (missing(maxSubjects)) maxSubjects=unique_check(sim.setup,'maxSubjects')
  if (missing(dropout_rate)) dropout_rate=unique_check(sim.setup,'dropout_rate')
  if (missing(sf_rate)) sf_rate=unique_check(sim.setup,'sf_rate')
  if (missing(accrual_rate)) accrual_rate=unique_check(sim.setup,'accrual_rate')
  if (missing(window_screen)) window_screen=unique_check(sim.setup,'window_screen')
  if (missing(window_DLT)) window_DLT=unique_check(sim.setup,'window_DLT')
  if (missing(meanDLTtime)) meanDLTtime=unique_check(sim.setup,'meanDLTtime')
  if (missing(lagPlanned)) lagPlanned=unique_check(sim.setup,'lagPlanned')
  if (missing(lagDropout)) lagDropout=unique_check(sim.setup,'lagDropout')
  if (missing(k.earlystop)) k.earlystop=unique_check(sim.setup,'k.earlystop')
  if (missing(n.earlystop)) n.earlystop=unique_check(sim.setup,'n.earlystop')
  
  targetRate=unique(sim.setup$targetRate)
  p.saf=unique(sim.setup$p.saf)
  p.tox=unique(sim.setup$p.tox)
  dose=1:length(sim.setup$trueRate[[1]])
  
  rows_to_plot = which(sim.setup$cohortSize==cohortSize &
                         sim.setup$min_cohortEnrolled==min_cohortEnrolled &
                         sim.setup$min_cohortEval==min_cohortEval &
                         sim.setup$maxSubjects==maxSubjects &
                         sim.setup$dropout_rate==dropout_rate &
                         sim.setup$sf_rate==sf_rate &
                         sim.setup$accrual_rate==accrual_rate &
                         sim.setup$window_screen==window_screen &
                         sim.setup$window_DLT==window_DLT &
                         sim.setup$meanDLTtime==meanDLTtime &
                         sim.setup$lagPlanned==lagPlanned &
                         sim.setup$lagDropout==lagDropout &
                         sim.setup$k.earlystop==k.earlystop &
                         sim.setup$n.earlystop==n.earlystop)
  
  new_object=lapply(object,function(df) df[rows_to_plot])
  sim.setup=new_object$sim.setup%>%bind_rows()
  
  # true rates for all the toxicity scenarios included in the boinsim object
  trueRate=unique(sim.setup$trueRate)
  trueRate=t(as.data.frame(trueRate))
  rownames(trueRate)<-NULL
  colnames(trueRate)<-paste('Dose',dose)
  trueRate=as.data.frame(trueRate)
  
  df_trueRate=trueRate%>%mutate(Total=NA)%>%gather(dose,n)
  
  # compile operating characteristics for table generation
  if (dim(trueRate)[1]>1){ # more than one toxicity scenario is used
    df_SlcPerc = new_object$SlcPerc_dose%>%bind_rows()%>%select(-`All excluded`)%>%mutate(Total=rowSums(.))%>%gather(dose,n)
    df_DLTs = new_object$DLTs%>%bind_rows()%>%gather(dose,n)%>%mutate(perc=100*n/new_object$SubjectsCompleted%>%bind_rows()%>%gather(dose,perc)%>%.$perc)
    df_subjectsRect = new_object$Subjects%>%bind_rows()%>%gather(dose,n)%>%mutate(perc=new_object$Subjects_perc%>%bind_rows()%>%gather(dose,perc)%>%.$perc)
    df_subjectsEnrl = new_object$SubjectsEnrolled%>%bind_rows()%>%gather(dose,n)%>%mutate(perc=new_object$SubjectsEnrolled_perc%>%bind_rows()%>%gather(dose,perc)%>%.$perc)
    df_subjectsEval = new_object$SubjectsCompleted%>%bind_rows()%>%gather(dose,n)%>%mutate(perc=new_object$SubjectsCompleted_perc%>%bind_rows()%>%gather(dose,perc)%>%.$perc)
  }
  else { # only one toxicity scenario is included 
    # one-hot encoding for dose selected at each simulation run to transform it into the same data frame format as other operating characteristics
    tmp = data.frame(dose=new_object$DoseSelected)
    tmp$dose<-factor(tmp$dose,levels=c(0,dose),labels = c('All excluded',paste('Dose',dose)))
    dvar=dummyVars(~dose,data=tmp,fullRank=FALSE,levelsOnly=T)
    df_SlcPerc=as.data.frame(predict(dvar, newdata = tmp))%>%select(-`All excluded`)%>%
      mutate(Total=rowSums(.))%>%
      mutate(across(.fns=as.character))%>%
      gather(dose,n)%>%mutate(across(.cols=n,.fns=as.numeric))%>%
      group_by(dose)%>%summarise(n=sum(n)/dim(sim.setup)[1]*100,.groups='drop')
    
    df_DLTs = new_object$DLTs%>%bind_rows()%>%mutate(Total=rowSums(.))%>%gather(dose,n)%>%mutate(perc=100*n/new_object$SubjectsCompleted%>%bind_rows()%>%mutate(Total=rowSums(.))%>%gather(dose,perc)%>%.$perc)
    df_subjectsRect = new_object$Subjects%>%bind_rows()%>%mutate(Total=rowSums(.))%>%gather(dose,n)%>%mutate(perc=new_object$Subjects_perc%>%bind_rows()%>%mutate(Total=rowSums(.))%>%gather(dose,perc)%>%.$perc)
    df_subjectsEnrl = new_object$SubjectsEnrolled%>%bind_rows()%>%mutate(Total=rowSums(.))%>%gather(dose,n)%>%mutate(perc=new_object$SubjectsEnrolled_perc%>%bind_rows()%>%mutate(Total=rowSums(.))%>%gather(dose,perc)%>%.$perc)
    df_subjectsEval = new_object$SubjectsCompleted%>%bind_rows()%>%mutate(Total=rowSums(.))%>%gather(dose,n)%>%mutate(perc=new_object$SubjectsCompleted_perc%>%bind_rows()%>%mutate(Total=rowSums(.))%>%gather(dose,perc)%>%.$perc)
  }
  
  # trial duration summary data frame
  df_duration = data.frame(`Trial duration` = new_object$DurationTrial, `Avg.cohort duration`=new_object$avg_DurationCohort)
  
  # Table 1
  df_list=list(df_trueRate,df_SlcPerc,df_DLTs,df_subjectsRect,df_subjectsEnrl,df_subjectsEval)
  names=c('True Toxicity Rate','% MTD selection','Number (rate) of DLTs ','Number (%) of subjects recruited','Number (%) of subjects enrolled','Number (%) of subjects evaluable')
  ## helper function: calculate summary statistics
  summarystats<-function(df){
    if ('perc'%in% colnames(df)){
      summary<-df%>%group_by(dose)%>%
        summarise(n_mean=signif(mean(n,na.rm=TRUE),3),n_sd=signif(sd(n,na.rm=TRUE),3),
                  n_median=signif(median(n,na.rm=TRUE),3),n_Q1=signif(quantile(n,0.25,na.rm=TRUE),3),n_Q3=signif(quantile(n,0.75,na.rm=TRUE),3),
                  perc_mean=signif(mean(perc,na.rm=TRUE),3),perc_sd=signif(sd(perc,na.rm=TRUE),3),
                  perc_median=signif(median(perc,na.rm=TRUE),3),perc_Q1=signif(quantile(perc,0.25,na.rm=TRUE),3),perc_Q3=signif(quantile(perc,0.75,na.rm=TRUE),3),.groups='drop')
      }
    else {
      summary<-df%>%group_by(dose)%>%
        summarise(n_mean=signif(mean(n,na.rm=TRUE),3),n_sd=signif(sd(n,na.rm=TRUE),3),
                  n_median=signif(median(n,na.rm=TRUE),3),n_Q1=signif(quantile(n,0.25,na.rm=TRUE),3),n_Q3=signif(quantile(n,0.75,na.rm=TRUE),3),.groups='drop')
    }
    return(summary)
  }
  
  df_summary_list=lapply(df_list,summarystats)
  names(df_summary_list)=names
  
  ## helper function: format table output
  tbl_format<-function(df){
    if (dim(df%>%select(starts_with('perc')))[2]!=0){
      tbl.out = df %>% mutate(Mean=paste(n_mean,'(',paste0(perc_mean,'%'),')'),
                              SD=paste(n_sd,'(',paste0(perc_sd,'%'),')'),
                              Median=paste(n_median,'(',paste0(perc_median,'%'),')'),
                              IQR=paste(paste0('[',n_Q1,',',n_Q3,']'),
                                        '(',
                                        paste0('[',paste0(perc_Q1,'%'),',',paste0(perc_Q3,'%'),']'),
                                        ')'))%>%
        select(dose,Mean,SD,Median,IQR)%>%
        pivot_longer(Mean:IQR,names_to = 'Statistic',values_to = 'value')%>%
        pivot_wider(names_from = 'dose',values_from = 'value')%>%
        na.omit()
    }
    else{
      tbl.out = df %>% mutate(Mean=n_mean,
                              SD=n_sd,
                              Median=n_median,
                              IQR=paste('[',n_Q1,',',n_Q3,']'))%>%
        mutate(across(.fns=as.character))%>%
        select(dose,Mean,SD,Median,IQR)%>%
        pivot_longer(Mean:IQR,names_to = 'Statistic',values_to = 'value')%>%
        pivot_wider(names_from = 'dose',values_from = 'value')%>%
        na.omit()
    }
    if (!'SD' %in% tbl.out$Statistic) tbl.out=tbl.out%>%filter(Statistic=='Mean')
    return(tbl.out)
  }
  
  tbl.out_list=lapply(df_summary_list,tbl_format)
  
  tbl.out_list = lapply(seq_along(tbl.out_list), 
         function(i,list,names){
           new_df=list[[i]]%>%mutate(`Operating characteristics`=names[i])%>%relocate(`Operating characteristics`)
           if (grepl('subjects',names[i])){
             new_df = new_df%>%rowwise()%>%mutate(Total=strsplit(Total,split=" ")[[1]][1])
           }
           return(new_df)
         },list=tbl.out_list,names=names(tbl.out_list))
  
  
  tbl.out1 = bind_rows(tbl.out_list)
  table1 = flextable::flextable(tbl.out1)%>%
    merge_v(j=1)%>%
    align(align='center',part='all')%>%
    padding(padding=8,part='all')%>%
    set_table_properties(width=1,layout='autofit')
  hline_index = (1:dim(tbl.out1)[1])[tbl.out1$`Operating characteristics`!=lead(tbl.out1$`Operating characteristics`)]
  hline_index=hline_index[!is.na(hline_index)]
  table1<-hline(table1,i=hline_index)
  table1<-bold(table1, bold = TRUE, part = "header")
  table1<-footnote(table1,i=1,j=dim(tbl.out1)[2],
                  value=as_paragraph('When total % MTD selection is not 100%, 100-Total is the % of "All excluded".'),
                  ref_symbols = "1",part='header')
  table1<-set_caption(table1,caption='Operating characteristics summary table')
  
  #Table 2: duration
  df_duration_summary = df_duration%>%gather(measure,value)%>%group_by(measure,.drop=F)%>%
    summarise(Mean=signif(mean(value,na.rm=TRUE),3),SD=signif(sd(value,na.rm=TRUE),3),Median=signif(median(value,na.rm=TRUE),3),Q1=signif(quantile(value,0.25,na.rm=TRUE),3),Q3=signif(quantile(value,0.75,na.rm=TRUE),3),.groups='drop')
  tbl.out2 = df_duration_summary%>% mutate(IQR=paste('[',Q1,',',Q3,']'))%>%
    mutate(across(.fns=as.character))%>%
    select(measure,Mean,SD,Median,IQR)%>%
    pivot_longer(Mean:IQR,names_to = 'Statistic',values_to = 'value')%>%
    pivot_wider(names_from = 'measure',values_from = 'value')%>%
    na.omit()
  table2 = flextable::flextable(tbl.out2)%>%
    align(align='center',part='all')%>%
    padding(padding=8,part='all')%>%
    set_table_properties(width=1,layout='autofit')
  table2<-bold(table2, bold = TRUE, part = "header")
  table2<-set_caption(table2,caption='Trial duration (days) summary table')
  
  if (tableOutput){
    print(table1)
    print(table2)
  }
  else{
    outputlist = df_summary_list
    outputlist$Duration = df_duration_summary
    return(outputlist)
  }
}
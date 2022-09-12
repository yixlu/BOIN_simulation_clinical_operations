#' Generate operating characteristics plots.
#' 
#' Generate plots for true toxicity rates, dose-wise or interval-wise number of DLTs, number of recruited/enrolled/evaluable subjects, and MTD selection percentage.
#' 
#' @usage plotOC_boinsim(object,type,which=c(1,2,5,6),
#'               cohortSize,min_cohortEnrolled,min_cohortEval,maxSubjects,
#'               dropout_rate,sf_rate,accrual_rate,
#'               window_screen,window_DLT,meanDLTtime,lagPlanned,lagDropout,
#'               k.earlystop,n.earlystop)
#' @param object an object output by run_boinsim or gridsearch_boinsim. 
#' @param type 'dose' or 'interval'. Whether to plot by groups of doses or intervals (below, within, or above the target interval).
#' @param which a subset of the numbers 1:6, by default 1,2,5,6, referring to:
#'              1. true toxicity rate at each dose level (a line plot if only one toxicity scenario is include, a jitterplot otherwise);
#'              2. # of DLTs
#'              3. # of subjects recruited
#'              4. # of subjects enrolled
#'              5. # of subjects evaluable
#'              6. MTD selection percentage
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
#' 
#' @details - If \code{object} is an output by run_boinsim, \code{rawOutput} has to be set to TRUE in run_boinsim. The plots will be a summary of all simulation trials.
#'            There is no need to specify parameters other than \code{type} and \code{which}.
#'          - If \code{object} is an output by gridsearch_boinsim, the plots will be a summary across all toxicity scenarios, but only for one combination of other parameters as specified by the rest of the arguments.
#'            For the arguments that only take one value in \code{object} (not subject to gridsearch), there is no need to specify here again. 
#' 
#' @import dplyr tidyr ggplot2 scales caret
#' @export

plotOC_boinsim<- function(object,type,which=c(1,2,5,6),
                   cohortSize,min_cohortEnrolled,min_cohortEval,maxSubjects,
                   dropout_rate,sf_rate,accrual_rate,
                   window_screen,window_DLT,meanDLTtime,lagPlanned,lagDropout,
                   k.earlystop,n.earlystop){
  
  sim.setup=object$sim.setup%>%bind_rows()
  ## helper function
  unique_check<-function(sim.setup,parameter){
    tmp=sim.setup%>%pull(parameter)
    tmp=unique(tmp)
    if (length(tmp)!=1){stop('Only one value is allowed for each parameter.')}
    return(tmp)
  }
  # filter dataset with specified parameters
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
  
  df_trueRate=trueRate%>%
    gather(dose,trueRate)%>%
    rowwise()%>%
    mutate(TargetInterval=(trueRate>p.saf)+(trueRate>p.tox)) %>%
    mutate(TargetInterval=factor(TargetInterval,levels=c(0,1,2),labels = c('Below target interval','Within target interval','Above target interval') )) %>%
    ungroup()
  
  # true toxicity rate plot
  if (dim(trueRate)[1]==1){ # plot as line plot if only one toxicity scenario is included
    p1<-ggplot2::ggplot(df_trueRate,aes(x=dose,y=trueRate,group=1))+
      geom_line()+
      geom_point()+
      geom_hline(aes(yintercept=targetRate,linetype='target'),color='black')+
      geom_hline(aes(yintercept =p.saf,linetype='p.saf'),color='blue')+
      geom_hline(aes(yintercept =p.tox,linetype='p.tox'),color='red')+
      scale_linetype_manual(name = "Target Interval", values = c(3,3,2),guide = guide_legend(override.aes = list(color = c("black","blue", "red"))))+
      ylim(0,1)+
      theme_classic()+
      theme(legend.position = 'bottom')
  }
  else {
    p1<-ggplot2::ggplot(df_trueRate,aes(x=dose,y=trueRate))+
      geom_violin()+
      geom_jitter(aes(color=TargetInterval),shape=16,position=position_jitter(0.2))+
      scale_color_manual(breaks = c('Below target interval','Within target interval','Above target interval'),
                         values=scales::muted(c("blue", "green", "red"),60,95))+
      geom_hline(aes(yintercept=targetRate,linetype='target'),color='black')+
      geom_hline(aes(yintercept =p.saf,linetype='p.saf'),color='blue')+
      geom_hline(aes(yintercept =p.tox,linetype='p.tox'),color='red')+
      scale_linetype_manual(name = "Target Interval", values = c(3,3,2),guide = guide_legend(override.aes = list(color = c("black","blue", "red"))))+
      ylim(0,1)+
      theme_classic()+
      theme(legend.position = 'bottom',legend.title = element_blank())
  }
  
  # compile operating characteristics for plot generation based on choice of 'type'
  if (type=='dose'){
    df_compile=data.frame(level=rep(paste('Dose',dose),each=dim(sim.setup)[1]),id=1:(length(dose)*dim(sim.setup)[1]))
    if (dim(trueRate)[1]>1){ # more than one toxicity scenario is used
      df_SlcPerc = new_object$SlcPerc_dose%>%bind_rows()%>%select(-`All excluded`)
      df_DLTs = new_object$DLTs%>%bind_rows()%>%select(-Total)
      df_subjectsRect = new_object$Subjects%>%bind_rows()%>%select(-Total)
      df_subjectsEnrl = new_object$SubjectsEnrolled%>%bind_rows()%>%select(-Total)
      df_subjectsEval = new_object$SubjectsCompleted%>%bind_rows()%>%select(-Total)
    }
    else { # only one toxicity scenario is included 
      # one-hot encoding for dose selected at each simulation run to transform it into the same data frame format as other operating characteristics
      tmp = data.frame(dose=new_object$DoseSelected)
      tmp$dose<-factor(tmp$dose,levels=c(0,dose),labels = c('All excluded',paste('Dose',dose)))
      dvar=dummyVars(~dose,data=tmp,fullRank=FALSE,levelsOnly=T)
      df_SlcPerc=as.data.frame(predict(dvar, newdata = tmp))%>%select(-`All excluded`)
      
      df_DLTs = new_object$DLTs%>%bind_rows()
      df_subjectsRect = new_object$Subjects%>%bind_rows()
      df_subjectsEnrl = new_object$SubjectsEnrolled%>%bind_rows()
      df_subjectsEval = new_object$SubjectsCompleted%>%bind_rows()
    }
  }
  else if (type=='interval'){
    df_compile=data.frame(level=rep(c('Below target interval','Within target interval','Above target interval'),each=dim(sim.setup)[1]),id=1:(3*dim(sim.setup)[1]))
    if (dim(trueRate)[1]>1){
      df_SlcPerc = new_object$SlcPerc_interval%>%bind_rows()%>%select(-`All excluded`)
      df_DLTs = new_object$DLTs_interval%>%bind_rows()%>%select(-Total)
      df_subjectsRect = new_object$SubjectsInterval%>%bind_rows()%>%select(-Total)
      df_subjectsEnrl = new_object$SubjectsEnrolledInterval%>%bind_rows()%>%select(-Total)
      df_subjectsEval = new_object$SubjectsCompletedInterval%>%bind_rows()%>%select(-Total)
    }
    else {
      tmp = data.frame(dose=new_object$DoseSelected_interval)
      dvar=dummyVars(~dose,data=tmp,fullRank=FALSE,levelsOnly=T)
      df_SlcPerc=as.data.frame(predict(dvar, newdata = tmp))%>%select(-`All excluded`)
      
      df_DLTs = new_object$DLTs_interval%>%bind_rows()
      df_subjectsRect = new_object$SubjectsInterval%>%bind_rows()
      df_subjectsEnrl = new_object$SubjectsEnrolledInterval%>%bind_rows()
      df_subjectsEval = new_object$SubjectsCompletedInterval%>%bind_rows()
    }
  }
  
  df_list=list(df_SlcPerc,df_DLTs,df_subjectsRect,df_subjectsEnrl,df_subjectsEval)
  names=c('SlcPerc','DLTs','subjectsRect','subjectsEnrl','subjectsEval')
  for (i in 1:length(df_list)){
    name=names[i]
    df_compile=df_list[[i]]%>%
      gather(level,!!name)%>%
      mutate(id=1:dim(df_compile)[1])%>%
      left_join(df_compile,by=c('id','level'))
  }
  df_compile<-df_compile%>%select(-id)
  
  
  if (type=='dose'){
    # # of DLTs
    p2<-ggplot2::ggplot(df_compile,aes(x=level,y=DLTs))+
      geom_boxplot()+
      ylab('# of DLTs')+
      xlab('Dose level')+
      scale_y_continuous(expand = expansion(c(0,.2),0))+
      theme_bw()+
      theme(panel.grid.major.x = element_blank(),panel.grid.minor = element_blank(),
            axis.text= element_text(size=14), axis.title = element_text(size=14,face='bold'))
    # # of subjectsRect
    p3<-ggplot2::ggplot(df_compile,aes(x=level,y=subjectsRect))+
      geom_boxplot()+
      ylab('# of subjects recruited')+
      xlab('Dose level')+
      scale_y_continuous(expand = expansion(c(0,.2),0))+
      theme_bw()+
      theme(panel.grid.major.x = element_blank(),panel.grid.minor = element_blank(),
            axis.text= element_text(size=14), axis.title = element_text(size=14,face='bold'))
    # # of subjectsEnrl
    p4<-ggplot2::ggplot(df_compile,aes(x=level,y=subjectsEnrl))+
      geom_boxplot()+
      ylab('# of subjects enrolled')+
      xlab('Dose level')+
      scale_y_continuous(expand = expansion(c(0,.2),0))+
      theme_bw()+
      theme(panel.grid.major.x = element_blank(),panel.grid.minor = element_blank(),
            axis.text= element_text(size=14), axis.title = element_text(size=14,face='bold'))
    # # of subjectsEval
    p5<-ggplot2::ggplot(df_compile,aes(x=level,y=subjectsEval))+
      geom_boxplot()+
      ylab('# of subjects evaluable')+
      xlab('Dose level')+
      scale_y_continuous(expand = expansion(c(0,.2),0))+
      theme_bw()+
      theme(panel.grid.major.x = element_blank(),panel.grid.minor = element_blank(),
            axis.text= element_text(size=14), axis.title = element_text(size=14,face='bold'))
    # MTD selection percentage
    if (dim(trueRate)[1]==1){ # barplot if only one toxicity scenario
      SlcPerc=df_compile%>%group_by(level,.drop=F)%>%summarize(slcp=sum(SlcPerc)/dim(sim.setup)[1]*100,.groups = 'drop')
      p6<-ggplot2::ggplot(SlcPerc,aes(x=level,y=slcp))+
        geom_bar(stat='identity')+
        geom_text(data=SlcPerc,aes(x=level,y=slcp,label=paste0(slcp,'%')),vjust=-0.2)+
        ylab('MTD selection percentage')+
        xlab('Dose level')+
        scale_y_continuous(expand = expansion(c(0,.2),0))+
        theme_bw()+
        theme(panel.grid.major.x = element_blank(),panel.grid.minor = element_blank(),
              axis.text= element_text(size=14), axis.title = element_text(size=14,face='bold'))
    }
    else { # boxplot if more than one toxicity scenario
      p6<-ggplot2::ggplot(df_compile,aes(x=level,y=SlcPerc))+
        geom_boxplot()+
        ylab('MTD selection percentage')+
        xlab('Dose level')+
        scale_y_continuous(expand = expansion(c(0,.2),0))+
        theme_bw()+
        theme(panel.grid.major.x = element_blank(),panel.grid.minor = element_blank(),
              axis.text= element_text(size=14), axis.title = element_text(size=14,face='bold'))
    }
  }
  else if (type=='interval'){
    df_compile$level<-factor(df_compile$level,levels=c("Below target interval", "Within target interval", "Above target interval"))
    # # of DLTs
    p2<-ggplot2::ggplot(df_compile,aes(x=level,y=DLTs,fill=level))+
      geom_boxplot()+
      scale_fill_manual(breaks = c("Below target interval", "Within target interval", "Above target interval"),
                        values=scales::muted(c("blue", "green", "red"),60,95))+
      ylab('# of DLTs')+
      xlab('Location on dose-toxicity curve')+
      scale_y_continuous(expand = expansion(c(0,.2),0))+
      theme_bw()+
      theme(panel.grid.major.x = element_blank(),panel.grid.minor = element_blank(),
            axis.text= element_text(size=14), axis.title = element_text(size=14,face='bold'),
            legend.position = 'bottom',legend.title = element_blank())
    
    # # of subjectsRect
    p3<-ggplot2::ggplot(df_compile,aes(x=level,y=subjectsRect,fill=level))+
      geom_boxplot()+
      scale_fill_manual(breaks = c("Below target interval", "Within target interval", "Above target interval"),
                        values=scales::muted(c("blue", "green", "red"),60,95))+
      ylab('# of subjects recruited')+
      xlab('Location on dose-toxicity curve')+
      scale_y_continuous(expand = expansion(c(0,.2),0))+
      theme_bw()+
      theme(panel.grid.major.x = element_blank(),panel.grid.minor = element_blank(),
            axis.text= element_text(size=14), axis.title = element_text(size=14,face='bold'),
            legend.position = 'bottom',legend.title = element_blank())
    
    # # of subjectsEnrl
    p4<-ggplot2::ggplot(df_compile,aes(x=level,y=subjectsEnrl,fill=level))+
      geom_boxplot()+
      scale_fill_manual(breaks = c("Below target interval", "Within target interval", "Above target interval"),
                        values=scales::muted(c("blue", "green", "red"),60,95))+
      ylab('# of subjects enrolled')+
      xlab('Location on dose-toxicity curve')+
      scale_y_continuous(expand = expansion(c(0,.2),0))+
      theme_bw()+
      theme(panel.grid.major.x = element_blank(),panel.grid.minor = element_blank(),
            axis.text= element_text(size=14), axis.title = element_text(size=14,face='bold'),
            legend.position = 'bottom',legend.title = element_blank())
    
    # # of subjectsEval
    p5<-ggplot2::ggplot(df_compile,aes(x=level,y=subjectsEval,fill=level))+
      geom_boxplot()+
      scale_fill_manual(breaks = c("Below target interval", "Within target interval", "Above target interval"),
                        values=scales::muted(c("blue", "green", "red"),60,95))+
      ylab('# of subjects evaluable')+
      xlab('Location on dose-toxicity curve')+
      scale_y_continuous(expand = expansion(c(0,.2),0))+
      theme_bw()+
      theme(panel.grid.major.x = element_blank(),panel.grid.minor = element_blank(),
            axis.text= element_text(size=14), axis.title = element_text(size=14,face='bold'),
            legend.position = 'bottom',legend.title = element_blank())
    # MTD selection percentage
    if (dim(trueRate)[1]==1){ # barplot if only one toxicity scenario
      SlcPerc=df_compile%>%group_by(level,.drop=F)%>%summarize(slcp=sum(SlcPerc)/dim(sim.setup)[1]*100,.groups = 'drop')
      p6<-ggplot2::ggplot(SlcPerc,aes(x=level,y=slcp,fill=level))+
        geom_bar(stat='identity')+
        scale_fill_manual(breaks = c("Below target interval", "Within target interval", "Above target interval"),
                          values=scales::muted(c("blue", "green", "red"),60,95))+
        geom_text(data=SlcPerc,aes(x=level,y=slcp,label=paste0(slcp,'%')),vjust=-0.2)+
        ylab('MTD selection percentage')+
        xlab('Location on dose-toxicity curve')+
        scale_y_continuous(expand = expansion(c(0,.2),0))+
        theme_bw()+
        theme(panel.grid.major.x = element_blank(),panel.grid.minor = element_blank(),
              axis.text= element_text(size=14), axis.title = element_text(size=14,face='bold'),
              legend.position = 'bottom',legend.title = element_blank())
    }
    else { # boxplot if more than one toxicity scenario
      p6<-ggplot2::ggplot(df_compile,aes(x=level,y=SlcPerc,fill=level))+
        geom_boxplot()+
        ylab('MTD selection percentage')+
        xlab('Location on dose-toxicity curve')+
        scale_fill_manual(breaks = c("Below target interval", "Within target interval", "Above target interval"),
                          values=scales::muted(c("blue", "green", "red"),60,95))+
        scale_y_continuous(expand = expansion(c(0,.2),0))+
        theme_bw()+
        theme(panel.grid.major.x = element_blank(),panel.grid.minor = element_blank(),
              axis.text= element_text(size=14), axis.title = element_text(size=14,face='bold'),
              legend.position = 'bottom',legend.title = element_blank())
    }
  }
  plots_list=list(p1,p2,p3,p4,p5,p6)
  plots_list[which]
}
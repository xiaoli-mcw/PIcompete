#' Sample-Weighted Prevalence-Incidence Mixture Models for Competing Incident Events (as of 03/03/2020)
#' 
#' This package fits competing risks models to failure time (or survival time) data for two competing events. 
#' Failure time for one event of the competing events can be prevalent left-censored, interval-censored or a mixture 
#' of truly incident disease and missed prevalent disease when disease ascertainment is not always conducted at baseline, 
#' while failure time for the other event is only interval-censored. Baseline is set to be time 0. 
#' General transformation,G(x)=(1+r*x)/r if r>0; =x if r=0, is used for a subdistribution hazard function multiplied by an exponential effect of a linear combination of risk factors  
#' for flexible incidence models. Logistic regression models are used for prevalence. The IPW log-likelihood approach, which uses the inverse of sample inclusion probabilities,
#' is employed to account for different sampling fractions across strata. 
#' 
#' @import fdrtool nloptr
#' 
#' @param Data Data used to fit the model containing columns for each term in p.model, i.model1 and i.model2 expressions.
#'  For stratified random sampling designs, columns denoted samp.wgt and strata are expected indicating the sampling weights and sampling strata.
#'  population="super" option, an additional column denoted strata.frac is expected indicating the fraction of the population
#'  that consists of each strata.  For example, if in the target population there are three strata that occurs with proportions 0.2, 0.4, and 0.6,
#'  then strata.frac will take values of 0.2, 0.4 or 0.6.
#' @param p.model The prevalence model for event 1 to be fitted, specified using an expression of the form \emph{C~model}.
#'  Elements in the expression are as followed:
#'  \itemize{
#'  \item c - Numeric variable indicating whether the event was prevalent at time zero, taking values of 1="Yes", 0="No", -999="Unknown";
#'  \item model - Linear predictor consisting of a series of terms separated by \emph{+} operators.
#'  }
#' @param i.model1 The incidence model for event 1 to be fitted using an expression of the form \emph{C+L1+R1~model1}
#'  \itemize{
#'  \item C - Numeric variable indicating whether the event was prevalent at time zero, taking values of 1="Yes", 0="No", -999="Unknown";
#'  \item L1 - Numeric starting time of the interval in which event 1 occurred, with -999 denoting known prevalent events;
#'  \item R1 - Ending time of the interval in which event 1 occurred, with -999 and Inf denoting known prevalent event 1 and right-censoring, respectively;
#'  \item model1 - Linear predictor consisting of a series of terms separated by \emph{+} operators.
#'  }
#' @param i.model2 The incidence model for event 2 to be fitted, specified using an expression of the form \emph{L2+R2~model2}
#' \itemize{
#'  \item L2 - Numeric starting time of the interval in which event 1 occurred, with -999 denoting known prevalent event 1;
#'  \item R2 - Ending time of the interval in which event 1 occurred, with -999 and Inf denoting known prevalent event 1 and right-censoring, respectively;
#'  \item model2 - Linear predictor consisting of a series of terms separated by \emph{+} operators.
#'  }
#' @param trans.r1  The parameter "r" for the transformation function for event 1, G(x)=log(1+rx)/r for r>0;G(x)=x for r=0 (default),which indicates proportional hazards model for the subdistribution hazard function.
#' @param trans.r2  The parameter "r" for the transformation function for event 2. Default to 0. 
#' @param n.beta is The number of regressors expressed in the p.model plus 1 (for intercept). If p.model is "C~1", n.beta=1.
#' @param n.gamma1  The number of regressors expressed in the i.model1. If i.model1 is "C+L1+R1~1", n.gamma1=0.
#' @param n.gamma2 The number of regressors expressed in the i.model2. If i.model2 is "L2+R2~1", n.gamma2=0.
#' @param reg.initials The initial values for regression coefficients in the order of (p.model, i.model1, i.model2). The number of components for reg.initials is n.beta+n.gamma1+n.gamma2. Default to be NULL.
#' @param convergence.criteria The criterion for the convergence of the iterated algorithm. Default to 0.001
#' @param iteration.limit The maximum number allowed for the iteration of the algorithm. Default to 250.
#' @param time.interval time.interval determines how finner finite time points are evenly divided, at which subdistribution hazard functions are estimated. Default to 0.1.
#' @param time.list a vector of finite time points at which subdistribution hazard functions are estimated. Default to NULL. For example, when an irregular spaced time points are of interest, time.list=c(1,3,8,10).
#' @param population options="super" and "finite" include variation due to super-population sampling and finite sampling from the super-population and variation due to finite sampling from a finite population, respectively. Default to "super".
#' @param anal.var analytical variance estimation is provided when anal.var=TRUE and the inverse information matrix exists. Default to TRUE.
#'
#' @return The output is a list of class ipw.pi.competing.risks, which contains the following elements.
#'  \itemize{
#'  \item data.summary: A data frame containing the following:
#'    Included subjects - number of observations with complete data;
#'    Known prevalent event 1 - the number of events known to be prevalent at time zero;
#'    Incident event 1 - the number of event times for event 1 occuring in the interval (L1>=0,R1<Inf] and C=0;
#'    Incident event 2 - the number of event times for event 2 occuring in the interval (L2>=0,R2<Inf] and C=0;
#'    Left censored event 1 - the number of event times known to occur by R1<Inf, but can also have been prevalent at time zero, that is C=-999;
#'    Right censoring - the number of observations right-censored with event time occurring in the interval (L1>0,R1=Inf) or (L2>0,R2=Inf) with C=0;
#'    Missing prevalent+right censoring - the number of observations with intervals (0,Inf) for both events 1 and 2 and C=-999.
#'    Maximum follow-up time for event 1
#'    Maximum follow-up time for event 2
#'    Maximum right censoring time for event 1
#'    Maximum right censoring time for event 2
#'  \item reg.coef: A data frame summarizing parameter values, standard errors, and 95 percent confidence intervals.
#'  \item reg.covariance: The analytical asymptotic covariance matrix exists for the regression coefficient estimates
#'  \item prevalence: A vector of prevalence estimates given covariates specified in p.model 
#'  \item subdist.hazard1: A data frame includes two columns, time and estimated subdistribution hazard for event 1.
#'  \item subdist.hazard2: A data frame includes two columns, time and estimated subdistribution hazard for event 2.
#'  \item subdist.hazard.fn1: A function returns estimated subdistribution hazard for event 1 when time points are input.
#'  \item subdist.hazard.fn2: A function returns estimated subdistribution hazard for event 2 when time points are input.
#'  \item convergence: Convergence statistics
#'  \item run.time.mins: The elapsed computation time in minutes.
#'  \item loglike: Sample-weighted log-likelikelihood at the estimated parameters
#'  \item trans.r: A vector of the specified parameters for the transformation functions used for events 1 and 2
#'  \item models: A vector of the specified models, p.model, i.model1 and i.model2
#'  }
#'
#' @author Noorie Hyun, \email{nhyun@mcw.edu}, Xiao Li \email{xiaoli@mcw.edu} 
#'
#' @references
#' \itemize{
#'  \item Hyun N, Katki HA, Graubard BI.
#'        Sample-Weighted Semiparametric Estimation of Cause-Specific Cumulative Risk and Incidence Using Left or Interval-Censored Data from Electronic Health Records. Statistics in Medicine 2020;
#'        under the 2nd review.
#' }
#'
#' @export
#'


ipw.pi.competing<-function(Data,p.model,i.model1,i.model2,trans.r1=0,trans.r2=0,n.beta=1, n.gamma1=0,n.gamma2=0,
                                 reg.initials=NULL,convergence.criteria=0.001,iteration.limit=250,
                                 time.interval=0.1,time.list=NULL,population="super",anal.var=TRUE,...){
  #library('fdrtool') is required for Lambda estimates
  #library('nloptr') is required for regression parameter estimates with constraints
  #population option "super" or "finite"
  ########## Transformation functions #############################
  trans<-function(r,x){
    if(r==0){
      output<-x
    }else if(r>0){
      output<-log(1+r*x)/r
    }
    return(output)
  }
  
  dG<-function(r,x){
    if(r==0){
      output<-1
    }else if(r>0){
      output<-1/(1+r*x)
    }
    return(output)
  }
  
  inv.trans<-function(r,x){
    if(r==0){
      output<-x
    }else if(r>0){
      output<- (exp(r*x)-1)/r
    }
    return(output)
  }
  ######################### semiparametric weighted log-like ########################################################
  wobs.semipara.like<-function(na.drop=FALSE,trans.r1,trans.r2,para,n.beta,n.gamma1,n.gamma2,design.mat,xmat1,xmat2,sdata3){
    #sdata3 should include the variables of LambdaL and LambddaR
    
    n.para<-n.beta+n.gamma1+n.gamma2
    mm<-n.beta+1
    mm2<-n.beta+n.gamma1
    mm3<-n.beta+n.gamma1+1
    
    betta<-para[1:n.beta]
    if(n.gamma1>0&n.gamma2>0){
      gam1<-para[mm:mm2]
      gam2<-para[mm3:n.para]
      gamma.risk1<-xmat1%*%gam1
      gamma.risk2<-xmat2%*%gam2
    } else if(n.gamma1==0&n.gamma2>0){
      gam1<-0
      gam2<-para[mm3:n.para]
      gamma.risk1<-0
      gamma.risk2<-xmat2%*%gam2
    } else if(n.gamma1>0&n.gamma2==0){
      gam1<-para[mm:mm2]
      gam2<-0
      gamma.risk1<-xmat1%*%gam1
      gamma.risk2<-0
    }else if(n.gamma1==0&n.gamma2==0){
      gam1<-0
      gam2<-0
      gamma.risk1<-0
      gamma.risk2<-0
    }
    
    if(n.beta>1){
      beta.risk<-design.mat%*%betta
    }else if(n.beta==1){
      beta.risk<-design.mat*betta
    }
    
    
    LambdaL1<-sdata3$LambdaL1
    LambdaR1<-sdata3$LambdaR1
    LambdaL2<-sdata3$LambdaL2
    LambdaR2<-sdata3$LambdaR2
    event<-sdata3$event
    C<-sdata3$C
    samp.wgt<-sdata3$samp.wgt
    
    # 3 cases
    #C=-999 and event=1 or -999
    #C=-999 and event=2/3 ;C=0 and event=1/2/3
    #C=1
    
    #event==1/-999
    gp1<-which(C==0&event==1)
    gp2<-which(C==0&event==2)
    gp3<-which(C==0&event==3)
    gp4<-which(C==-999&event==1)
    gp5<-which(C==-999&event==3)
    
    
    prev<-exp(beta.risk)/(1+exp(beta.risk))
    
    
    lp<-log(prev)*samp.wgt
    lcomp<-log(1-prev)*samp.wgt
    
    
    e.risk1<-as.numeric(exp(gamma.risk1))
    e.risk2<-as.numeric(exp(gamma.risk2))
    
    ######### Transformation of Lambda(t)*e.risk #################
    
    G.L1<-trans(trans.r1,LambdaL1*e.risk1)
    G.R1<-trans(trans.r1,LambdaR1*e.risk1)
    G.L2<-trans(trans.r2,LambdaL2*e.risk2)
    G.R2<-trans(trans.r2,LambdaR2*e.risk2)
    
    
    cond1<-exp(-G.L1)-exp(-G.R1)
    cond2<-exp(-G.L2)-exp(-G.R2)
    cond3<-exp(-G.L1)+exp(-G.L2)
    sgr1<-which(cond1>0)
    sgr2<-which(cond2>0)
    sgr3<-which(cond3>1)
    
    
    B.1<-log(    (exp(-G.L1)-exp(-G.R1))[intersect(gp1,sgr1)] )*samp.wgt[intersect(gp1,sgr1)]
    B.2<-log(    (exp(-G.L2)-exp(-G.R2))[intersect(gp2,sgr2)])*samp.wgt[intersect(gp2,sgr2)]

    if(na.drop==TRUE){
      B.3<-log(    (exp(-G.L1)+exp(-G.L2)-1 )[intersect(gp3,sgr3)])*samp.wgt[intersect(gp3,sgr3)]
    }else if(na.drop==FALSE){
      B.3<-log(    (exp(-G.L1)+exp(-G.L2)-1 )[gp3])*samp.wgt[gp3]
    }
    #B.3[is.na(B.3)|B.3==-Inf]<-0
    
    B.4<-log( (prev+(1-prev)*(1-exp(-G.R1)))[gp4])*samp.wgt[gp4]
    B.5<-log( (prev+(1-prev)*exp(-G.L2))[gp5])*samp.wgt[gp5]
    
    output<-sum(lp[C==1],lcomp[C==0],B.1,B.2,B.3,B.4,B.5)
    
    return(output)
  }  #pseudo-log-likelihood
  
  
  wsemi.opt <- function(para){
    -wobs.semipara.like(na.drop=FALSE,trans.r1,trans.r2,para,n.beta,n.gamma1,n.gamma2,design.mat,xmat1,xmat2,sdata3)
  }
  wsemi.opt.subgr <- function(para){
    -wobs.semipara.like(na.drop=TRUE,trans.r1,trans.r2,para,n.beta,n.gamma1,n.gamma2,design.mat,xmat1,xmat2,sdata3)
  }
  
  
  score.f<-function(trans.r1,trans.r2,para,n.beta,n.gamma1,n.gamma2,design.mat,xmat1,xmat2,sdata3){
    
    n.samp<-nrow(sdata3)
    n.para<-n.beta+n.gamma1+n.gamma2
    mm<-n.beta+1
    mm2<-n.beta+n.gamma1
    mm3<-n.beta+n.gamma1+1
    
    betta<-para[1:n.beta]
   
    if(n.beta>1){
      beta.risk<-design.mat%*%betta
    }else if(n.beta==1){
      beta.risk<-design.mat*betta
    }
    
    if(n.gamma1>0){
      gam1<-para[mm:mm2]
      gamma.risk1<-xmat1%*%gam1
    } else if(n.gamma1==0){
      gam1<-0
      gamma.risk1<-0
    }
    
    if(n.gamma2>0){
      gam2<-para[mm3:n.para]
      gamma.risk2<-xmat2%*%gam2
    } else if(n.gamma2==0){
      gam2<-0
      gamma.risk2<-0
    }
    
    
    
    LambdaL1<-sdata3$LambdaL1
    LambdaR1<-sdata3$LambdaR1
    LambdaL2<-sdata3$LambdaL2
    LambdaR2<-sdata3$LambdaR2
    event<-sdata3$event
    C<-sdata3$C
    samp.wgt<-sdata3$samp.wgt
    
    # 3 cases
    #C=-999 and event=1 or -999
    #C=-999 and event=2/3 ;C=0 and event=1/2/3
    #C=1
    
    prev<-exp(beta.risk)/(1+exp(beta.risk))
    comp.prev<-1/(1+exp(beta.risk))
    
    lp<-log(prev)
    lcomp<-log(1-prev)
    
    
    e.risk1<-as.numeric(exp(gamma.risk1))
    e.risk2<-as.numeric(exp(gamma.risk2))
    e.beta.risk<-as.numeric(exp(beta.risk))
    
    
    gp1<-which(C==0&event==1)
    gp2<-which(C==0&event==2)
    gp3<-which(C==0&event==3)
    gp4<-which(C==-999&event==1)
    gp5<-which(C==-999&event==3)
    
    L1.e1<-LambdaL1*e.risk1
    R1.e1<-LambdaR1*e.risk1
    L2.e2<-LambdaL2*e.risk2
    R2.e2<-LambdaR2*e.risk2
    
    G.L1<-trans(trans.r1,L1.e1)
    G.R1<-trans(trans.r1,R1.e1)
    G.L2<-trans(trans.r2,L2.e2)
    G.R2<-trans(trans.r2,R2.e2)
    
    
    S2.L<-exp(- G.L2)
    S2.R<-exp(-G.R2)
    S1.L<-exp(-G.L1)
    S1.R<-exp(-G.R1)
    
    #gradient of beta
    gr.beta<-matrix(rep(n.samp*n.beta),n.samp,n.beta)
    beta.cons<-rep(0,n.samp)
    beta.cons[C==1]<-comp.prev[C==1] #C==1
    beta.cons[C==0]<-(-prev)[C==0]
    beta.cons[gp5]<-(prev*(1-S2.L)/(e.beta.risk+S2.L))[gp5]
    beta.cons[gp4]<-(prev*S1.R/(e.beta.risk+1-S1.R))[gp4]
    gr.beta<-beta.cons*design.mat*samp.wgt
    gr.beta.nowt<-beta.cons*design.mat
    
    #gradient of gamma1
    if(n.gamma1>0){
      gamma1.cons<-rep(0,n.samp)
      gamma1.cons[gp1]<-((S1.R*dG(trans.r1,R1.e1)*R1.e1-S1.L*dG(trans.r1,L1.e1)*L1.e1)/(S1.L-S1.R) )[gp1]
      gamma1.cons[gp3]<-(-S1.L*dG(trans.r1,L1.e1)*L1.e1/(S1.L+S2.L-1) )[gp3]
      gamma1.cons[gp4]<-(S1.R*dG(trans.r1,R1.e1)*R1.e1/(e.beta.risk+1-S1.R) )[gp4]
      
      gr.gamma1<-xmat1*gamma1.cons*samp.wgt
      gr.gamma1.nowt<-xmat1*gamma1.cons
    }
    if(n.gamma2>0){
      gamma2.cons<-rep(0,n.samp)
      gamma2.cons[gp2]<-((S2.R*dG(trans.r2,R2.e2)*R2.e2-S2.L*dG(trans.r2,L2.e2)*L2.e2)/(S2.L-S2.R) )[gp2]
      gamma2.cons[gp3]<-((-S2.L*dG(trans.r2,L2.e2)*L2.e2)/(S1.L+S2.L-1) )[gp3]
      gamma2.cons[gp5]<-(-S2.L*dG(trans.r2,L2.e2)*L2.e2/(e.beta.risk+S2.L) )[gp5]
      
      gr.gamma2<-xmat2*gamma2.cons*samp.wgt
      gr.gamma2.nowt<-xmat2*gamma2.cons
    }
     
    if(n.gamma1>0&n.gamma2>0){
      output<-cbind(gr.beta,gr.gamma1,gr.gamma2)
      output.nowt<-cbind(gr.beta.nowt,gr.gamma1.nowt,gr.gamma2.nowt)
    }else if(n.gamma1>0&n.gamma2==0){
        output<-cbind(gr.beta,gr.gamma1)
        output.nowt<-cbind(gr.beta.nowt,gr.gamma1.nowt)
    }else if(n.gamma1==0&n.gamma2>0){
      output<-cbind(gr.beta,gr.gamma2)
      output.nowt<-cbind(gr.beta.nowt,gr.gamma2.nowt)
    }else if(n.gamma1==0&n.gamma2==0){
        output<-gr.beta
        output.nowt<-gr.beta.nowt
    }  
    
    output22<-list()
    output[is.na(output)]<-0
    output.nowt[is.na(output.nowt)]<-0
    
    output22$score.vec<-output
    output22$score.vec.nowt<-output.nowt
    output22$score<-colSums(output,na.rm=TRUE)
    return(output22)
  } #score.f
  
  #In the simulation setting prep1_5 with no covariate, this initial values work well.
  sub.haz.ini2<-function(x){
   hazard.ini<-0.15*(1-exp(-0.15*x))
    return(hazard.ini)
  }
  
  
  #In the simulation setting prep2_1; this initial values work well, in particular, w.r.t gamma2.1 estimates and Lambda 2 estimate.
  sub.haz.ini1<-function(x){
    hazard.ini<-0.25*(1-exp(-0.25*x))
    return(hazard.ini)
  }
  
  
  grad.f<-function(para){
    gradient<-score.f(trans.r1,trans.r2,para,n.beta,n.gamma1,n.gamma2,design.mat,xmat1,xmat2,sdata3)
    return(-gradient$score)
  }
  
  
  cov.mat.ph1<-function(delta.theta,samp.data){
    n.para<-ncol(delta.theta)
    cov.mat<-0*diag(n.para)
    
    aaa<-which(table(samp.data$strata)!=1)
    strata.list<-as.numeric(names(aaa))
    
    for (k in strata.list){
      average<-apply(delta.theta[samp.data$strata==k,],2,mean)
      delta.theta.gp<-delta.theta[samp.data$strata==k,]
      n.gp<-nrow(delta.theta.gp)
      ave.mat<-matrix(average,n.gp,n.para,1)
      c.delta.theta.gp<-delta.theta.gp-ave.mat
      cov.mat<-cov.mat+n.gp/(n.gp-1)*(t(c.delta.theta.gp)%*%c.delta.theta.gp)
    }
    return(cov.mat)
  }#the end of the function cov.mat.ph1
  
  
  
  cov.mat.ph2<-function(score.nowt,samp.data){
    n.para<-ncol(score.nowt)
    cov.mat<-0*diag(n.para)
    
    aaa<-which(table(samp.data$strata)!=1)
    strata.list<-as.numeric(names(aaa))
    
    for (k in strata.list){
      
      score.by.strata<-score.nowt[samp.data$strata==k,]
      n.gp<-nrow(score.by.strata)
      strata.fraction<-as.numeric(samp.data$strata.frac[samp.data$strata==k])[1]
      samp.data$s.frc<-1/samp.data$samp.wgt
      samp.fraction<-as.numeric(samp.data$s.frc[samp.data$strata==k])[1]
      A<- (1/n.gp)*(t(score.by.strata)%*%score.by.strata)
      B<-colMeans(score.by.strata)
      B2<-B%*%t(B)
      cond.var<-A-B2
      cov.mat<-cov.mat+(strata.fraction*(1-samp.fraction)/samp.fraction)*(cond.var)
    }
    return(cov.mat)
  }#the end of the function cov.mat.ph2
  
  
  Lambda.initials<-function(time.list,Lambda,missing.upper,Data,event,time.1,time.2,cap.Lambda){
    #Data may not include LambdaL and LambdaR
    if(event==1){
      
      #Lambda1<-time.list1*0.015
      Lambda.data<-data.frame(time=time.list,Lambda=Lambda)
      
      Lambda.index<-rbind(Lambda.data,rep(0,2),c(Inf,cap.Lambda),rep(missing.upper,2))
      sdata<-Data
      
      sdata2<-merge(sdata,Lambda.index,by.x="L1",by.y="time",all.x=TRUE)
      colnames(sdata2)[ncol(sdata2)]<-"LambdaL1"
      sdata3<-merge(sdata2,Lambda.index,by.x="R1",by.y="time",all.x=TRUE)
      colnames(sdata3)[ncol(sdata3)]<-"LambdaR1"
      
      sdata3$LambdaL1[sdata3$L1<time.1]<-0
      sdata3$LambdaL1[sdata3$L1>time.2&sdata3$L1<missing.upper]<-max(sdata3$LambdaL1[sdata3$LambdaL1<missing.upper],na.rm=TRUE)
      sdata3$LambdaR1[sdata3$R1>time.2&sdata3$R1<missing.upper]<-max(sdata3$LambdaL1[sdata3$LambdaL1<missing.upper],na.rm=TRUE)
      
    }else if(event==2){
      Lambda.data<-data.frame(time=time.list,Lambda=Lambda)
      
      Lambda.index<-rbind(Lambda.data,rep(0,2),c(Inf,cap.Lambda),rep(missing.upper,2))
      sdata<-Data
      
      sdata2<-merge(sdata,Lambda.index,by.x="L2",by.y="time",all.x=TRUE)
      colnames(sdata2)[ncol(sdata2)]<-"LambdaL2"
      sdata3<-merge(sdata2,Lambda.index,by.x="R2",by.y="time",all.x=TRUE)
      colnames(sdata3)[ncol(sdata3)]<-"LambdaR2"
      
      sdata3$LambdaL2[sdata3$L2<time.1]<-0
      sdata3$LambdaL2[sdata3$L2>time.2&sdata3$L2<missing.upper]<-max(sdata3$LambdaL2[sdata3$LambdaL2<missing.upper],na.rm=TRUE)
      sdata3$LambdaR2[sdata3$R2>time.2&sdata3$R2<missing.upper]<-max(sdata3$LambdaL2[sdata3$LambdaL2<missing.upper],na.rm=TRUE)
    }
    output<-list()
    output$sdata3<-sdata3
    output$Lambda.data<-Lambda.data
    output$Lambda.index<-Lambda.index
    return(output)
  }
  
  
  
  L1update<-function(trans.r1,trans.r2,Data,tt.list,Lambda.data,missing.upper){
    
    #Data was updated with new regression parameters estimation (Data<-sdata3)
    #Data should include LambdaL1, LambdaR1, LambdaL2 and LambdaR2
    
    #sdata<-samp.data
    ssdata<-Data
    n.samp<-nrow(ssdata)
    
    
    ssdata$WL1<-rep(0,n.samp)
    ssdata$WL2<-rep(0,n.samp)
    ssdata$WL3<-rep(0,n.samp)
    ssdata$WL4<-rep(0,n.samp)
    samp.wgt<-ssdata$samp.wgt
    
    exc1<-which(abs(ssdata$LambdaR1-ssdata$LambdaL1)<1e-10)
    exc2<-which(is.na(ssdata$LambdaL1))
    exc3<-which(is.na(ssdata$LambdaR1))
    
    condition<-exp(-trans(trans.r1,ssdata$LambdaL1*ssdata$HR1) )+exp(-trans(trans.r2,ssdata$LambdaL2*ssdata$HR2))
    
    exc4<-which( condition<=1)
    
    
    subset1<-which(ssdata$C==0&ssdata$event==1&ssdata$L1>0)
    subset2<-which(ssdata$C==0&ssdata$event==1&ssdata$R1<missing.upper)
    subset3<-which(ssdata$C==0&ssdata$event==3&ssdata$L1>0)
    subset4<-which(ssdata$C==-999&ssdata$event==1&ssdata$R1<missing.upper)
    
    
    
    # weighted process
    WL1<- -exp(-trans(trans.r1,ssdata$LambdaL1*ssdata$HR1) )*ssdata$HR1*dG(trans.r1,ssdata$LambdaL1*ssdata$HR1)/(exp(-trans(trans.r1,ssdata$LambdaL1*ssdata$HR1) )-exp(-trans(trans.r1,ssdata$LambdaR1*ssdata$HR1) ))
    ssdata$WL1[subset1]<-WL1[subset1]*samp.wgt[subset1]
    ssdata$WL1[union(exc1,exc2)]<-0
    
    
    # weighted process
    WL2<- exp(-trans(trans.r1,ssdata$LambdaR1*ssdata$HR1) )*ssdata$HR1*dG(trans.r1,ssdata$LambdaR1*ssdata$HR1)/(exp(-trans(trans.r1,ssdata$LambdaL1*ssdata$HR1) )-exp(-trans(trans.r1,ssdata$LambdaR1*ssdata$HR1) ))
    ssdata$WL2[subset2]<-WL2[subset2]*samp.wgt[subset2]
    ssdata$WL2[union(exc1,exc3)]<-0
    
    # weighted process
    WL3<- -exp(-trans(trans.r1,ssdata$LambdaL1*ssdata$HR1) )*ssdata$HR1*dG(trans.r1,ssdata$LambdaL1*ssdata$HR1)/(exp(-trans(trans.r1,ssdata$LambdaL1*ssdata$HR1) )+exp(-trans(trans.r2,ssdata$LambdaL2*ssdata$HR2) )-1)
    ssdata$WL3[subset3]<-WL3[subset3]*samp.wgt[subset3]
    ssdata$WL3[union(exc4,exc2)]<-0
    
    # weighted process
    WL4<-exp(-trans(trans.r1,ssdata$LambdaR1*ssdata$HR1) )*ssdata$HR1*dG(trans.r1,ssdata$LambdaR1*ssdata$HR1)/(ssdata$expplus-exp(-trans(trans.r1,ssdata$LambdaR1*ssdata$HR1) ))
    ssdata$WL4[subset4]<-WL4[subset4]*samp.wgt[subset4]
    
    
    ssdata$WL1.sq<-ssdata$WL1^2
    ssdata$WL2.sq<-ssdata$WL2^2
    ssdata$WL3.sq<-ssdata$WL3^2
    ssdata$WL4.sq<-ssdata$WL4^2
    
    k<-0
    Q_L<-0
    W_L<-0
    G_L<-0
    TL<-nrow(Lambda.data)
    Q_Lambda<-dG_Lambda<-G_Lambda<-W_Lambda<-rep(NA,TL)
    
    for (ordered.t in Lambda.data$time){
      k<-k+1
      L.process<-as.numeric(ssdata$L1<=ordered.t)
      R.process<-as.numeric(ssdata$R1<=ordered.t)
      
      W_Lambda[k]<-sum(R.process*(ssdata$WL2+ssdata$WL4)+L.process*(ssdata$WL1+ssdata$WL3),na.rm=TRUE)
      G_Lambda[k]<-sum(R.process*(ssdata$WL2.sq+ssdata$WL4.sq)+L.process*(ssdata$WL1.sq+ssdata$WL3.sq),na.rm=TRUE)
      
      
      if(k==1){dG_Lambda[k]<-G_Lambda[k]
      } else if(k>1){
        dG_Lambda[k]<-G_Lambda[k]-G_Lambda[k-1]
      }
      Q_L<- Q_L+Lambda.data$Lambda[k]*dG_Lambda[k]
      Q_Lambda[k]<-W_Lambda[k]+Q_L
    }
    
    
    no.change<-which(dG_Lambda==0)
    if(length(no.change)>0){
      G_Lambda<-G_Lambda[-no.change]
      Q_Lambda<-Q_Lambda[-no.change]
      Lambda.data<-Lambda.data[-no.change,]
    }
    
    
    G_Lambda<-c(0,G_Lambda)
    Q_Lambda<-c(0,Q_Lambda)
    
    
    if(length(G_Lambda)!=length(unique(G_Lambda))){
      stop("The algorithm for estimating subdistribution hazard function for event 1 does not converge.")
    }
    
    GCM = gcmlcm(G_Lambda, Q_Lambda)
    
    time.knots<-which(G_Lambda[-1] %in% GCM$x.knots[-1])
    
    update.Lambda<-data.frame(time=Lambda.data$time[time.knots],Lambda=GCM$slope.knots)
    update.Lambda$Lambda[update.Lambda$Lambda<0]<-0
    
    sfun.1  <- stepfun(update.Lambda$time, c(update.Lambda$Lambda[1],update.Lambda$Lambda), f = 1) #CADLAG function
    new.Lambda<-data.frame(time=tt.list,Lambda=sfun.1(tt.list))
    
    
    output<-list()
    output$new.Lambda<-new.Lambda
    output$sfun<-sfun.1
    return(output)
  }#L1update
  
  
  
  L2update<-function(trans.r1,trans.r2,Data,tt.list,Lambda.data,missing.upper){
    
    #Data was updated with new regression parameters estimation
    #Data should include LambdaL1, LambdaR1, LambdaL2 and LambdaR2
    
    
    #sdata<-samp.data
    ssdata<-Data
    n.samp<-nrow(ssdata)
    
    ssdata$WL1<-rep(0,n.samp)
    ssdata$WL2<-rep(0,n.samp)
    ssdata$WL3<-rep(0,n.samp)
    ssdata$WL4<-rep(0,n.samp)
    samp.wgt<-ssdata$samp.wgt
    
    
    exc1<-which(abs(ssdata$LambdaR2-ssdata$LambdaL2)<1e-10)
    exc2<-which(is.na(ssdata$LambdaL2))
    exc3<-which(is.na(ssdata$LambdaR2))
    
    condition<-exp(-trans(trans.r1,ssdata$LambdaL1*ssdata$HR1) )+exp(-trans(trans.r2,ssdata$LambdaL2*ssdata$HR2) )
    
    exc4<-which( condition<=1)
    
    
    subset1<-which(ssdata$C==0&ssdata$event==2&ssdata$L2>0)
    subset2<-which(ssdata$C==0&ssdata$event==2&ssdata$R2<missing.upper)
    subset3<-which(ssdata$C==0&ssdata$event==3&ssdata$L2>0)
    subset4<-which(ssdata$C==-999&ssdata$event==3&ssdata$L2>0)
    
    
    # weighted process
    WL1<- -exp(-trans(trans.r2,ssdata$LambdaL2*ssdata$HR2))*ssdata$HR2*dG(trans.r2,ssdata$LambdaL2*ssdata$HR2)/(exp(-trans(trans.r2,ssdata$LambdaL2*ssdata$HR2) )-exp(-trans(trans.r2,ssdata$LambdaR2*ssdata$HR2) ))
    ssdata$WL1[subset1]<-WL1[subset1]*samp.wgt[subset1]
    ssdata$WL1[union(exc1,exc2)]<-0
    
    # weighted process
    WL2<- exp(-trans(trans.r2,ssdata$LambdaR2*ssdata$HR2) )*ssdata$HR2*dG(trans.r2,ssdata$LambdaR2*ssdata$HR2)/(exp(-trans(trans.r2,ssdata$LambdaL2*ssdata$HR2) )-exp(-trans(trans.r2,ssdata$LambdaR2*ssdata$HR2) ))
    ssdata$WL2[subset2]<-WL2[subset2]*samp.wgt[subset2]
    ssdata$WL2[union(exc1,exc3)]<-0
    
    # weighted process
    WL3<- -exp(-trans(trans.r2,ssdata$LambdaL2*ssdata$HR2) )*ssdata$HR2*dG(trans.r2,ssdata$LambdaL2*ssdata$HR2)/(exp(-trans(trans.r1,ssdata$LambdaL1*ssdata$HR1) )+exp(-trans(trans.r2,ssdata$LambdaL2*ssdata$HR2) )-1)
    ssdata$WL3[subset3]<-WL3[subset3]*samp.wgt[subset3]
    ssdata$WL3[union(exc4,exc1)]<-0
    
    # weighted process
    WL4<- -exp(-trans(trans.r2,ssdata$LambdaL2*ssdata$HR2) )*ssdata$HR2*dG(trans.r2,ssdata$LambdaL2*ssdata$HR2)/(rep(-1,n.samp)+ssdata$expplus+exp(-trans(trans.r2,ssdata$LambdaL2*ssdata$HR2) ))
    ssdata$WL4[subset4]<-WL4[subset4]*samp.wgt[subset4]
    
    
    ssdata$WL1.sq<-ssdata$WL1^2
    ssdata$WL2.sq<-ssdata$WL2^2
    ssdata$WL3.sq<-ssdata$WL3^2
    ssdata$WL4.sq<-ssdata$WL4^2
    
    k<-0
    Q_L<-0
    W_L<-0
    G_L<-0
    TL<-length(tt.list)
    Q_Lambda<-dG_Lambda<-G_Lambda<-W_Lambda<-rep(NA,TL)
    
    for (ordered.t in Lambda.data$time){
      k<-k+1
      L.process<-as.numeric(ssdata$L2<=ordered.t)
      R.process<-as.numeric(ssdata$R2<=ordered.t)
      
      W_Lambda[k]<-sum(R.process*(ssdata$WL2)+L.process*(ssdata$WL1+ssdata$WL3+ssdata$WL4),na.rm=TRUE)
      G_Lambda[k]<-sum(R.process*(ssdata$WL2.sq)+L.process*(ssdata$WL1.sq+ssdata$WL3.sq+ssdata$WL4.sq),na.rm=TRUE)
      
      if(k==1){dG_Lambda[k]<-G_Lambda[k]
      } else if(k>1){
        dG_Lambda[k]<-G_Lambda[k]-G_Lambda[k-1]
      }
      Q_L<- Q_L+Lambda.data$Lambda[k]*dG_Lambda[k]
      Q_Lambda[k]<-W_Lambda[k]+Q_L
    }
    
    
    no.change<-which(dG_Lambda==0)
    if(length(no.change)>0){
      G_Lambda<-G_Lambda[-no.change]
      Q_Lambda<-Q_Lambda[-no.change]
      Lambda.data<-Lambda.data[-no.change,]
    }
    
    
    G_Lambda<-c(0,G_Lambda)
    Q_Lambda<-c(0,Q_Lambda)
    
    
    if(length(G_Lambda)!=length(unique(G_Lambda))){
      stop("The algorithm for estimating subdistribution hazard function for event 2 does not converge.")
    }
    
    GCM = gcmlcm(G_Lambda, Q_Lambda)
    
    time.knots<-which(G_Lambda[-1] %in% GCM$x.knots[-1])
    
    update.Lambda<-data.frame(time=Lambda.data$time[time.knots],Lambda=GCM$slope.knots)
    update.Lambda$Lambda[update.Lambda$Lambda<0]<-0
    
    
    sfun.1  <- stepfun(update.Lambda$time, c(update.Lambda$Lambda[1],update.Lambda$Lambda), f = 1) #CADLAG function
    new.Lambda<-data.frame(time=tt.list,Lambda=sfun.1(tt.list))
    
    #return(new.Lambda)
    
    output<-list()
    output$new.Lambda<-new.Lambda
    output$sfun<-sfun.1
    
    return(output)
    
  }#L2update
  
  
  wPI.get.summary<- function(samp.data) {
    n.samp<-nrow(samp.data)
    
    prevalent.event1<- sum(as.numeric(samp.data$C==1&samp.data$event==1))
    incident.event1<-sum(as.numeric(samp.data$C==0&samp.data$event==1))
    incident.event2<-sum(as.numeric(samp.data$C==0&samp.data$event==2))
    right.censoring.C0<-sum(as.numeric(samp.data$C==0&samp.data$event==3))
    right.censoring.C999<-sum(as.numeric(samp.data$C==-999&samp.data$event==3))
    left.censoring.event1<-sum(as.numeric(samp.data$C==-999&samp.data$event==1))
    max.event1.time<-max(samp.data$R1[samp.data$R1<Inf])
    max.event2.time<-max(samp.data$R2[samp.data$R2<Inf])
    max.right.censoring.time1<-max(samp.data$L1[samp.data$L1>0&samp.data$R1==Inf])
    max.right.censoring.time2<-max(samp.data$L2[samp.data$L2>0&samp.data$R2==Inf])
    
    
    
    data.summary1<-data.frame(c("Included subjects","Known prevalent event1","Incident event1","Incident event2",
                                "Left censored event1","Right censoring","Missing prevalent+right censoring", "maximum FU event1",
                                "maximum FU event2","maximum right censonring time event1","maximum right censonring time event2"),
                              c(n.samp,prevalent.event1,incident.event1,incident.event2,left.censoring.event1,
                                right.censoring.C0,right.censoring.C999,max.event1.time,max.event2.time,
                                max.right.censoring.time1,max.right.censoring.time2), stringsAsFactors=FALSE)
    
    names(data.summary1)<- c("label", "no.obs")
    
    return(data.summary1)
  }
  
  prev.fun<-function(beta){
    prev<-exp(beta)/(1+exp(beta))
    return(prev)
  }
  
  ################ Data Preparation ####################
  #change variable names for outcomes, L, R
  fml0<-as.formula(p.model)
  fml1<-as.formula(i.model1)
  fml2<-as.formula(i.model2)
  
  if(attr(terms(fml0),"response")!=1){
    stop("No response variable in model1")
  }
  
  if(attr(terms(fml1),"response")!=1){
    stop("No response variable in model1")
  }
  
  if(attr(terms(fml2),"response")!=1){
    stop("No response variable in model2")
  }
  if (!(c("event") %in% names(Data))) {
    stop("No event variable in the data")
  }
  
  allvars0<-all.vars(fml0)
  allvars1<-all.vars(fml1)
  allvars2<-all.vars(fml2)
  
  
  names(Data)[which(names(Data) %in% allvars1[1])]<-"C"
  names(Data)[which(names(Data) %in% allvars1[2])]<-"L1"
  names(Data)[which(names(Data) %in% allvars1[3])]<-"R1"
  
  names(Data)[which(names(Data) %in% allvars2[1])]<-"L2"
  names(Data)[which(names(Data) %in% allvars2[2])]<-"R2"
  
  allvars0[1]<-c("C")
  allvars1[1:3]<-c("C","L1","R1")
  allvars2[1:2]<-c("L2","R2")
  
  
  
  #change the response variable name in p.model formula
  p.model.split<-strsplit(p.model,split="~")
  p.model1<-paste("C~",p.model.split[[1]][2])
  fml0<-as.formula(p.model1)
  
  i.model1.split<-strsplit(i.model1,split="~")
  i.model11<-paste("C+L1+R1~",i.model1.split[[1]][2])
  fml1<-as.formula(i.model11)
  
  i.model2.split<-strsplit(i.model2,split="~")
  i.model22<-paste("L2+R2~",i.model2.split[[1]][2])
  fml2<-as.formula(i.model22)
  
  remove(p.model1,i.model11,i.model22)
  
  
  #change the response variable name in i.model formula
  if(c("samp.wgt") %in% names(Data)) {
    allvars<-c(allvars0,allvars1,allvars2,"samp.wgt")
    if(sum(Data$samp.wgt<0)>0){
      stop("\"samp.weight\" should be positive.")
    }
  } else if(!(c("samp.wgt") %in% names(Data))){
    allvars<-c(allvars0,allvars1,allvars2,"samp.wgt")
    Data$samp.wgt<-1
  }
  
  if(c("strata") %in% names(Data)) {
    allvars<-c(allvars,"strata")
  } else if(!(c("strata") %in% names(Data))){
    allvars<-c(allvars,"strata")
    Data$strata<-1
  }
  if(population=="super"){
    if(c("strata.frac") %in% names(Data)){
      allvars<-c(allvars,"strata.frac") 
    }
  }
  
  allvars<-c(allvars,"event")
  allvars<-unique(allvars)
  
  samp.data2<-Data[,allvars]
  samp.data<-samp.data2[complete.cases(samp.data2),]  #excluding cases with NA
  
  ###############################################################################################
  #time.points2<-time points to be set for cumulative hazard estimate
  #fixed.initials initial parameters for beta and gamma just in case when logistic model doesn't converge.
  
  n.samp<-nrow(samp.data)
  n.regpara<-n.beta+n.gamma1+n.gamma2
  n.gamma<-n.gamma1+n.gamma2
  
  mm<-n.beta+1
  mm2<-n.beta+n.gamma1
  mm3<-n.beta+n.gamma1+1
  
  iter.limit<-400                          #iteration for cumulative hazard function
  
  ########### Initial parameters  ################################################
  
  if(is.null(reg.initials)){
    
    if(n.gamma1==0&n.gamma2==0){
      if (n.beta==1){
        initials1<--3
        initials2<--3
      }else if(n.beta>1){
        initials1<-c(-3,rep(0.1,n.beta-1))
        
        ini.data<-samp.data
        ini.data<-ini.data[ini.data$C==1|ini.data$C==0,]
        beta.ini<-glm(fml0, data = ini.data, family = "binomial")
        initials2<-beta.ini$coefficients
        remove(beta.ini,ini.data)
      }
    } else if(n.gamma1>0&n.gamma2>0) {
      if(n.beta==1) {
        initials1<-c(-3,rep(0.1,n.gamma))
        ini.data<-samp.data
        ini.data$event1<-as.numeric(samp.data$event==1)
        ini.data$event2<-as.numeric(samp.data$event==2)
        gam1.ini<-glm(paste0("event1~",i.model1.split[[1]][2],sep=""), data = ini.data, family = "binomial")
        gam2.ini<-glm(paste0("event2~",i.model2.split[[1]][2],sep=""), data = ini.data, family = "binomial")
        
        initials2<-c(-3,gam1.ini$coefficients[-1],gam2.ini$coefficients[-1])
        remove(ini.data,gam1.ini,gam2.ini)
        
      }else if(n.beta>1){
        initials1<-c(-3,rep(0.1,n.beta-1),rep(0.1,n.gamma))
        
        ini.data<-samp.data
        ini.data$event1<-as.numeric(samp.data$event==1)
        ini.data$event2<-as.numeric(samp.data$event==2)
        gam1.ini<-glm(paste0("event1~",i.model1.split[[1]][2],sep=""), data = ini.data, family = "binomial")
        gam2.ini<-glm(paste0("event2~",i.model2.split[[1]][2],sep=""), data = ini.data, family = "binomial")
        
        ini.data<-ini.data[ini.data$C==1|ini.data$C==0,]
        beta.ini<-glm(fml0, data = ini.data, family = "binomial")
        
        initials2<-c(beta.ini$coefficients,gam1.ini$coefficients[-1],gam2.ini$coefficients[-1])
        remove(ini.data,beta.ini,gam1.ini,gam2.ini)
      }
    } else if(n.gamma1>0&n.gamma2==0){
      if(n.beta==1) {
        initials1<-c(-3,rep(0.1,n.gamma1))
        ini.data<-samp.data
        ini.data$event1<-as.numeric(samp.data$event==1)
        gam1.ini<-glm(paste0("event1~",i.model1.split[[1]][2],sep=""), data = ini.data, family = "binomial")
        
        initials2<-c(-3,gam1.ini$coefficients[-1])
        remove(ini.data,gam1.ini)
        
      }else if(n.beta>1){
        initials1<-c(-3,rep(0.1,n.beta-1),rep(0.1,n.gamma1))
        
        ini.data<-samp.data
        ini.data$event1<-as.numeric(samp.data$event==1)
        gam1.ini<-glm(paste0("event1~",i.model1.split[[1]][2],sep=""), data = ini.data, family = "binomial")
        
        ini.data<-ini.data[ini.data$C==1|ini.data$C==0,]
        beta.ini<-glm(fml0, data = ini.data, family = "binomial")
        
        initials2<-c(beta.ini$coefficients,gam1.ini$coefficients[-1])
        remove(ini.data,beta.ini,gam1.ini)
      }
    } else if(n.gamma1==0&n.gamma2>0){
      if(n.beta==1) {
        initials1<-c(-3,rep(0.1,n.gamma2))
        ini.data<-samp.data
        ini.data$event2<-as.numeric(samp.data$event==2)
        gam2.ini<-glm(paste0("event2~",i.model2.split[[1]][2],sep=""), data = ini.data, family = "binomial")
        
        initials2<-c(-3,gam2.ini$coefficients[-1])
        remove(ini.data,gam2.ini)
        
      }else if(n.beta>1){
        initials1<-c(-3,rep(0.1,n.beta-1),rep(0.1,n.gamma2))
        
        ini.data<-samp.data
        
        ini.data$event2<-as.numeric(samp.data$event==2)
        gam2.ini<-glm(paste0("event2~",i.model2.split[[1]][2],sep=""), data = ini.data, family = "binomial")
        
        ini.data<-ini.data[ini.data$C==1|ini.data$C==0,]
        beta.ini<-glm(fml0, data = ini.data, family = "binomial")
        
        initials2<-c(beta.ini$coefficients,gam2.ini$coefficients[-1])
        remove(ini.data,beta.ini,gam2.ini)
    }
  } else {
    if(length(reg.initials)!=n.regpara){
      stop("The number of regression parameter initials are not matched with the model.")
    } else {
      initials2<-reg.initials
    }
  }
  }
  missing.upper1<-max(samp.data$R1[samp.data$R1<Inf],samp.data$L1)*1.5+2 #event type=1;the initial Lambda=time.list*1.5 and 2 is an arbitrary positive number
  missing.upper2<-max(samp.data$R2[samp.data$R2<Inf],samp.data$L2)*1.5+2 #event type=2;the initial Lambda=time.list*1.5 and 2 is an arbitrary positive number
  missing.upper<-max(missing.upper1,missing.upper2)
  
  tab.event<-table(samp.data$event[samp.data$C==0],samp.data$samp.wgt[samp.data$C==0])
  tab.event2<-t(tab.event)*as.numeric(colnames(tab.event))
  marginal.prop<-colSums(prop.table(tab.event2))
  
  
  cap.subdist1<-marginal.prop[1]+marginal.prop[3]*0.4
  cap.subdist2<-marginal.prop[2]+marginal.prop[3]*0.4
  
  cap.L1<--log(1-cap.subdist1)
  cap.L2<--log(1-cap.subdist2)
  
  
  samp.data$L1[samp.data$L1==-999]<-missing.upper
  samp.data$R1[samp.data$R1==-999]<-missing.upper
  samp.data$L2[samp.data$L2==-999]<-missing.upper
  samp.data$R2[samp.data$R2==-999]<-missing.upper
  
  samp.data$id<-seq(1,nrow(samp.data))
  
  
  ############## Lambda1 Initials #################################################
  
  
  Jn1<-sort(unique(c(samp.data$L1[samp.data$L1>0&samp.data$L1<missing.upper],
                     samp.data$R1[samp.data$R1<missing.upper])))
  Jn2<-sort(unique(c(samp.data$L2[samp.data$L2>0&samp.data$L2<missing.upper],
                     samp.data$R2[samp.data$R2<missing.upper])))
  
  #ordered statistics for time: min{X_i: T_i<=X_i} and max{X_i:X_i<=T_i}
  time1.1<-min(samp.data$R1[samp.data$R1<missing.upper])
  time1.2<-max(samp.data$L1[samp.data$L1>0&samp.data$L1<missing.upper])
  
  time2.1<-min(samp.data$R2[samp.data$R2<missing.upper])
  time2.2<-max(samp.data$L2[samp.data$L2>0&samp.data$L2<missing.upper])
  
  if(time1.2==0|time2.2==0){
    stop("NPMLE condition: the maximum left time point should be positive.")
  }
  
  if(is.null(time.interval)){
    time.interval<- 0.1
    time.points1<-seq(0,time1.2,by=time.interval)
    time.points2<-seq(0,time2.2,by=time.interval)
  }else {
    time.points1<-seq(0,time1.2,by=time.interval)
    time.points2<-seq(0,time2.2,by=time.interval)
  }
  
  
  
  time.list1<-Jn1[Jn1>=time1.1&Jn1<=time1.2]
  time.list2<-Jn2[Jn2>=time2.1&Jn2<=time2.2]
  
  time1.95<-quantile(time.list1,c(0.95))
  time2.95<-quantile(time.list2,c(0.95))
  
  
    Lamb.ini1<-Lambda.initials(time.list1,sub.haz.ini2(time.list1),missing.upper,samp.data,event=1,time1.1,time1.2,cap.L1)
    old.Lambda1<-Lamb.ini1$Lambda.data
    Lamb.ini2<-Lambda.initials(time.list2,sub.haz.ini2(time.list2),missing.upper,Data=Lamb.ini1$sdata3,event=2,time2.1,time2.2,cap.L2)
    old.Lambda2<-Lamb.ini2$Lambda.data
    sdata3<-Lamb.ini2$sdata3
  
  
  
  mf0 <- model.frame(formula=fml0, data=sdata3)
  design.mat <- model.matrix(attr(mf0, "terms"), data=mf0)
  
  
  mf <- model.frame(formula=fml1, data=sdata3)
  xmat1 <- model.matrix(attr(mf, "terms"), data=mf)
  
  mf2 <- model.frame(formula=fml2, data=sdata3)
  xmat2 <- model.matrix(terms(fml2), data=mf2)
  
  if(ncol(xmat2)>1){
    xmat2<-as.matrix(xmat2[,-1])
  } else if(ncol(xmat2)==1){
    xmat2<-0
  }
  if(ncol(xmat1)>1){
    xmat1<-as.matrix(xmat1[,-1])
  } else if(ncol(xmat1)==1){
    xmat1<-0
  }
  
  if(n.beta>1){
    expplus.est<-as.numeric(1+exp(design.mat%*%initials1[1:n.beta]))
  }else if(n.beta==1){
    expplus.est<-as.numeric(1+exp(design.mat*initials1[1]))
  }
  
  if(n.gamma1>0){
    cox.hr.est1<-as.numeric(exp(xmat1%*%initials1[mm:mm2]))
  } else {
    cox.hr.est1<-1
  }
  if(n.gamma2>0){
    cox.hr.est2<-as.numeric(exp(xmat2%*%initials1[mm3:n.regpara]))
  } else {
    cox.hr.est2<-1
  }
  
  sdata3$HR1<-cox.hr.est1
  sdata3$HR2<-cox.hr.est2
  sdata3$expplus<-expplus.est
  
  ########################### LOOPING ###########################################################
  old.wNR.est<-initials2
  iteration<-0
  keep.going<-1
  ptm<-proc.time()
  diff.estimates<-100
  
  
  
  while(keep.going==1){
    
    iteration<- iteration +1
    
    ######################   Cumulative Hazard Estimate   ##############################
    # update Lambda2
    iter<-0
    keep.looping<-1
    old.Lambda22<-old.Lambda2
    
    #L2OUT5<-list()
    while(keep.looping==1){
      iter<-iter+1
      
      L2out<-L2update(trans.r1,trans.r2,sdata3,time.list2,old.Lambda22,missing.upper)
     # L2OUT5[[iter]]<-L2out
      
      new.Lambda2<-L2out$new.Lambda
      
      diff.Lambda2<-max(abs(new.Lambda2$Lambda[time.list2<time2.95]-old.Lambda22$Lambda[time.list2<time2.95]))
      
      old.Lambda22<-new.Lambda2
      
      
      #update sdata3 with new.Lambda2
      remove.ind2<-which(names(sdata3) %in% c("LambdaL2","LambdaR2"))
      Lamb.ini2<-Lambda.initials(time.list2,new.Lambda2$Lambda,missing.upper,sdata3[,-remove.ind2],event=2,time2.1,time2.2,cap.L2)
      
      sdata3<-Lamb.ini2$sdata3
      
      
      if (diff.Lambda2<convergence.criteria|iter>iter.limit){
        keep.looping<-0
      }
    } #sdata3 is sorted by L2 and R2
    #For now, Lambda2 estimate is back to the estimate from the previous step:
    remove.ind2<-which(names(sdata3) %in% c("LambdaL2","LambdaR2"))
    Lamb.ini2<-Lambda.initials(time.list2,old.Lambda2$Lambda,missing.upper,sdata3[,-remove.ind2],event=2,time2.1,time2.2,cap.L2)
    sdata3<-Lamb.ini2$sdata3
    
    ## update 1 ##
    iter<-0
    keep.looping<-1
    old.Lambda11<-old.Lambda1
    
    while(keep.looping==1){
      iter<-iter+1
      L1out<-L1update(trans.r1,trans.r2,sdata3,time.list1,old.Lambda11,missing.upper)
      new.Lambda1<-L1out$new.Lambda
      
      
      
      diff.Lambda1<-max(abs(new.Lambda1$Lambda[time.list1<time1.95]-old.Lambda11$Lambda[time.list1<time1.95]))
      
      old.Lambda11<-new.Lambda1
      
      
      #update sdata3 with new.Lambda1
      remove.ind1<-which(names(sdata3) %in% c("LambdaL1","LambdaR1"))
      Lamb.ini1<-Lambda.initials(time.list1,new.Lambda1$Lambda,missing.upper,sdata3[,-remove.ind1],event=1,time1.1,time1.2,cap.L1)
      
      
      sdata3<-Lamb.ini1$sdata3
      
      
      if (diff.Lambda1<convergence.criteria|iter>iter.limit){
        keep.looping<-0
      }
    } #sdata3 is sorted by L1 and R1
    #Lambda2 estimate is now updated:
    remove.ind2<-which(names(sdata3) %in% c("LambdaL2","LambdaR2"))
    Lamb.ini2<-Lambda.initials(time.list2,new.Lambda2$Lambda,missing.upper,sdata3[,-remove.ind2],event=2,time2.1,time2.2,cap.L2)
    sdata3<-Lamb.ini2$sdata3

    
    
    mf0 <- model.frame(formula=fml0, data=sdata3)
    design.mat <- model.matrix(attr(mf0, "terms"), data=mf0)
    label0<-paste0("logit:",colnames(design.mat))
    
    mf <- model.frame(formula=fml1, data=sdata3)
    xmat1 <- model.matrix(attr(mf, "terms"), data=mf)
    
    mf2 <- model.frame(formula=fml2, data=sdata3)
    xmat2 <- model.matrix(terms(fml2), data=mf2)
    
   
    if(ncol(xmat2)>1){
      label2<-paste0("i.model2:",colnames(xmat2)[-1])
        xmat2<-as.matrix(xmat2[,-1])
    } else if(ncol(xmat2)==1){
      xmat2<-0
      label2<-NULL
    }
    if(ncol(xmat1)>1){
      label1<-paste0("i.model1:",colnames(xmat1)[-1])
      xmat1<-as.matrix(xmat1[,-1])
    } else if(ncol(xmat1)==1){
      xmat1<-0
      label1<-NULL
    }
    
    labels<-c(label0,label1,label2)
    
    r.cens.ind<-which(sdata3$C==0&sdata3$R1==Inf&sdata3$R2==Inf)
    
    if(n.gamma1>0){
      const.mat1<-as.matrix(xmat1[r.cens.ind,])
    } 
    if(n.gamma2>0) {
      const.mat2<-as.matrix(xmat2[r.cens.ind,])
    }
    
    
    
    ################# start  updating regression parameters
    
    
    const.ineq <- function(para) {
      grisk1<-as.numeric(const.mat1%*%para[mm:mm2])
      grisk2<-as.numeric(const.mat2%*%para[mm3:n.regpara])
      L.erisk1<-exp(grisk1)*sdata3$LambdaL1[r.cens.ind]
      L.erisk2<-exp(grisk2)*sdata3$LambdaL2[r.cens.ind]
      S1<-exp(-trans(trans.r1,L.erisk1))
      S2<-exp(-trans(trans.r2,L.erisk2))
      constr=S1+S2-1.001  #constr>=0
      return(constr)
    }
    
    const.ineq2 <- function(para) {
      grisk1<-as.numeric(const.mat1%*%para[mm:mm2])
      grisk2<-0
      L.erisk1<-exp(grisk1)*sdata3$LambdaL1[r.cens.ind]
      L.erisk2<-exp(grisk2)*sdata3$LambdaL2[r.cens.ind]
      S1<-exp(-trans(trans.r1,L.erisk1))
      S2<-exp(-trans(trans.r2,L.erisk2))
      constr=S1+S2-1.001  #constr>=0
      return(constr)
    }
    
    const.ineq3 <- function(para) {
      grisk1<-0
      grisk2<-as.numeric(const.mat2%*%para[mm3:n.regpara])
      L.erisk1<-exp(grisk1)*sdata3$LambdaL1[r.cens.ind]
      L.erisk2<-exp(grisk2)*sdata3$LambdaL2[r.cens.ind]
      S1<-exp(-trans(trans.r1,L.erisk1))
      S2<-exp(-trans(trans.r2,L.erisk2))
      constr=S1+S2-1.001  #constr>=0
      return(constr)
    }
    
    #resw<-solnl(X=old.wNR.est,objfun=wsemi.opt.subgr,confun=const.ineq)
    if(n.gamma1>0&n.gamma2>0){
      resw<-cobyla(old.wNR.est, fn=wsemi.opt.subgr, hin = const.ineq)
    }else if(n.gamma1>0&n.gamma2==0){
      resw<-cobyla(old.wNR.est, fn=wsemi.opt.subgr, hin = const.ineq2)
    }else if(n.gamma1==0&n.gamma2>0){
      resw<-cobyla(old.wNR.est, fn=wsemi.opt.subgr, hin = const.ineq3)
    }else if(n.gamma1==0&n.gamma2==0){
      resw<-optim(old.wNR.est,wsemi.opt,gr=grad.f,method="BFGS")
    }

    
############### insert the cut codes
    
    ################# end of updating regression parameters
    
    
    wNR.est<-resw$par
    #print(wNR.est)
    
    diff.reg.para<-max(abs(old.wNR.est-wNR.est))
    
    if(n.beta>1){
      expplus.est<-as.numeric(1+exp(design.mat%*%wNR.est[1:n.beta]))
    }else if(n.beta==1){
      expplus.est<-as.numeric(1+exp(design.mat*wNR.est[1]))
    }
    
    
    if(n.gamma1>0){
      cox.hr.est1<-as.numeric(exp(xmat1%*%wNR.est[mm:mm2]))
    } else {
      cox.hr.est1<-1
    }
    if(n.gamma2>0){
      cox.hr.est2<-as.numeric(exp(xmat2%*%wNR.est[mm3:n.regpara]))
    } else {
      cox.hr.est2<-1
    }
    
    sdata3$HR1<-cox.hr.est1
    sdata3$HR2<-cox.hr.est2
    sdata3$expplus<-expplus.est
    
    
    
    diff.estimates<-max(diff.reg.para,diff.Lambda1,diff.Lambda2)
    diff.reg.est<-diff.reg.para
    diff.Lambda.est1<-diff.Lambda1
    diff.Lambda.est2<-diff.Lambda2
    
    old.wNR.est<-wNR.est
    old.Lambda1<-new.Lambda1
    old.Lambda2<-new.Lambda2
    
    if(diff.estimates<convergence.criteria|iteration>iteration.limit){
      keep.going<-0
      run.time<-((proc.time()-ptm)[3])/60
    }
    
    
  } #while(keep.going==1){
  
  
  
  par.est<-wNR.est
  prev.est<-prev.fun(par.est)
  
  if(is.null(time.list)){
    Lambda.data1<-new.Lambda1
    Lambda.data2<-new.Lambda2
  }else{
    Lambda.data1<-data.frame(time=time.list,Lambda=L1out$sfun(time.list))
    Lambda.data2<-data.frame(time=time.list,Lambda=L2out$sfun(time.list))
  }
  
  wloglike<-wobs.semipara.like(na.drop=FALSE,trans.r1,trans.r2,par.est,n.beta,n.gamma1,n.gamma2,design.mat,xmat1,xmat2,sdata3)
   
  
  ############ Variance estimation for regression parameter estimates ###########
        if(anal.var==TRUE){
          score.reg<-score.f(trans.r1,trans.r2,par.est,n.beta,n.gamma1,n.gamma2,design.mat,xmat1,xmat2,sdata3)
          score.nowt<-score.reg$score.vec.nowt
          score.wt<-score.reg$score.vec
          
          info.mat<- (t(score.wt)%*%score.nowt)
          N<-(sum(sdata3$samp.wgt))
          
          
          inv.info2<-try(solve(info.mat,tol=1e-300),silent=TRUE)
          
          if('try-error' %in% class(inv.info2)){
            cat("inverse-matrix problem in info.mat\n") 
          }else {
            inv.info<-inv.info2 
            
            if(population=="super"){
              inv.info.mat<- inv.info*N #because of subgroup analysis, sum of sampling weights replaces N
              
              cov.mat.phase1<-inv.info.mat/N
              
              cov.mat.by.strata<-cov.mat.ph2(score.nowt,sdata3)
              
              cov.mat.phase2<-(inv.info.mat%*%cov.mat.by.strata%*%inv.info.mat)/N
              
              cov.mat<-cov.mat.phase1+cov.mat.phase2
              theta.se<-sqrt(diag(cov.mat))  
            } else if(population=="finite"){
              delta.s<-score.reg$score.vec
              
              cov.mat.theta<-cov.mat.ph1(delta.s,samp.data)
              
              cov.mat<-inv.info%*%cov.mat.theta%*%inv.info
              theta.se<-sqrt(diag(cov.mat))
            }
            
            para.95UL<-as.numeric(par.est+1.96*theta.se)
            para.95LL<-as.numeric(par.est-1.96*theta.se)
            
            
          } #if('try-error' %in% class(inv.info2)){
        }#if anal.var==TRUE
        
        
        
        
  
  output<-list()
  output$data.summary<-wPI.get.summary(samp.data)
  if(exists('theta.se')){
    output$reg.coef<-data.frame(cbind(label=labels,coef=par.est,se=theta.se,ll95=para.95LL,ul95=para.95UL) )
    rownames(output$reg.coef)<-NULL
    output$reg.covariance<-cov.mat
  }else {
    output$reg.coef<-data.frame(cbind(label=labels,coef=par.est))
    rownames(output$reg.coef)<-NULL
    output$reg.covariance<-NA
  }
  
  output$prevalence<-prev.est
  output$subdist.hazard1<-Lambda.data1
  output$subdist.hazard2<-Lambda.data2
  output$subdist.hazard.fn1<-L1out$sfun
  output$subdist.hazard.fn2<-L2out$sfun
  output$convergence<-c(iteration=iteration,converg.overall=diff.estimates, convg.diff.reg=diff.reg.est,convg.diff.L1=diff.Lambda.est1,convg.diff.L2=diff.Lambda.est2)
  output$run.time.mins<-run.time
  output$loglike<-wloglike
  output$trans.r<-c(r1=trans.r1,r2=trans.r2)
  output$models<-c(p.model=p.model,i.model1=i.model1,i.model2=i.model2)
    
  #Used optim algorithm in each iteration
  
  
  
  
  return(output)
}  # the end of the function compting.lc.semipara


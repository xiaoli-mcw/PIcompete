#' To generate the empirical distribution for estimated sub-distribution hazard functions and regression coefficients, 
#' we adopted the two-phase weighted bootstrap procedure (Hyun et al, reference). We fit ipw.pi.competing function from main.R to 
#' the bootstrap weighted data.The most arguments except for three arguments, n.boot, parallel, n.core are necessary for fitting "ipw.pi.competing" function.
#' 
#' @import fdrtool nloptr parallel doParallel
#' 
#' @param n.boot the number of replications for bootstrap resampling 
#' @param parallel if parallel=TRUE, n.core should be specified. Default to FALSE.
#' @n.core the nummber of cores for parallel computing
#'
#' @return The output is a list of class ipw.bootstrap.ph2, which contains the following elements.
#'  \itemize{
#'  \item reg.coef: estimated regression coefficients
#'  \item time finite time points at which subdistribution hazards are estimated
#'  \item subdist.hazard1 estimated subdistribution hazard function for event 1
#'  \item subdist.hazard2 estimated subdistribution hazard function for event 2
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
#'
ipw.bootstrap.ph2<-function(n.boot,Data,p.model,i.model1,i.model2,trans.r1=0,trans.r2=0,n.beta=1, n.gamma1=0,n.gamma2=0,
                                 reg.initials=NULL,convergence.criteria=0.001,iteration.limit=250,
                                 time.interval,time.list=NULL,parallel=FALSE,n.core=NULL,...){
  
########## Functions #############################
  
  srs.bootstrap.data<-function(Data){
        
        if(!(c("strata") %in% names(Data))){
          Data$strata<-1
        }
        
        if(!(c("samp.wgt") %in% names(Data))){
          Data$samp.wgt<-1
        }
        
        strata.list<-as.numeric(names(table(Data$strata)))
        sample.size<-as.numeric(table(Data$strata))
        wgt.list<-sapply(strata.list,function(x) Data$samp.wgt[Data$strata==x][1],simplify=TRUE)
        frac.rate<-1/wgt.list
        strata.size<-sample.size*wgt.list
        n.samp<-nrow(Data)
        
        samp.data<-Data  
        
        perturb1<-rep(1,n.samp)
        perturb2<-rep(0,n.samp)
        perturb2[samp.data$samp.wgt==1]<-1
        perturb.var<-frac.rate/(2-frac.rate)
        
        ### Perturbation at Phase 1############################
        pp<-lapply(1:length(strata.list),function(x)  rgamma(sample.size[x],shape=1/perturb.var[x],scale=perturb.var[x]))
        i<-0
        for(x in strata.list){
          i<-i+1
          index<-which(samp.data$strata==x)
          perturb1[index]<-as.numeric(pp[[i]])
        }
        
        ### Perturbation at Phase 2############################
        multiple.factor<-as.integer(wgt.list)
        remainder<-(wgt.list-multiple.factor)*sample.size
        
        stratum.add<-lapply(1:length(strata.list),function(x) sample(which(samp.data$strata==strata.list[x]),replace=FALSE,size=remainder[x]))
        stratum.ph<-lapply(1:length(strata.list),function(x) sample(c(rep(which(samp.data$strata==strata.list[x]),multiple.factor[x]),stratum.add[[x]]) ,size=sample.size[x],replace=FALSE))
        sampled.stratum<-lapply(1:length(strata.list), function(x) as.numeric(as.character(names(table(stratum.ph[[x]])))))
        
        for(x in 1:length(strata.list)){
          perturb2[sampled.stratum[[x]]]<-as.numeric(table(stratum.ph[[x]]))
        }
        
        samp.data$samp.wgt<-samp.data$samp.wgt*perturb1*perturb2
        
        return(samp.data)
  }
  
  bootstrap.result<-function(Data){
    resamp.data<-srs.bootstrap.data(Data)
    boot.out<-try(ipw.pi.competing(Data=resamp.data,p.model,i.model1,i.model2,trans.r1,trans.r2,n.beta, n.gamma1,n.gamma2,
                                   convergence.criteria,iteration.limit,
                                   time.interval,time.list,anal.var=FALSE),silent = TRUE)
    return(boot.out)
  }
  ##################################################
  if(parallel==TRUE&is.null(n.core)){
    stop('Specify the number of corse, "n.core" argument for parallel computing.')
  }
  
  
  mm1<-n.beta+1
  mm2<-n.beta+n.gamma1
  mm3<-n.beta+n.gamma1+1
  mm4<-n.beta+n.gamma1+n.gamma2
  
  b.iter<-0
  boot.reg<-boot.L1<-boot.L2<-res1<-NULL
  
      
    if(parallel==TRUE){
      library(parallel)
      library(doParallel)
      
      cl <- makeCluster(n.core)
      clusterExport(cl,c("fdrtool","nloptr","bootstrap.result","srs.bootstrap.data","ipw.pi.competing"))
      registerDoParallel(cl)
      res1 <- foreach(n = 1:n.boot*2, .combine = c, .packages = c("fdrtool","nloptr","bootstrap.result","srs.bootstrap.data","ipw.pi.competing"), .errorhandling=c('pass')) %dopar% 
        when(b.iter==0){b.out<-bootstrap.result(Data))
              if(is.matrix((summary(res1)))){
               success.list<-which(summary(res1)[,1]!=0)
               n.success<-length(success.list)
               if(n.success>=n.boot) {
                 b.iter<-1
               }
              }    
      }#when
      stopCluster(cl) 
      
      
       summary.reg.coef<-sapply(X=success.list,FUN=function(X) res1[[X]]$reg.para$coef)
       regpara.summary<-data.frame(cbind(t(summary.reg.coef[1:n.beta,]),
             t(summary.reg.coef[mm1:mm2,]),
             t(summary.reg.coef[mm3:mm4,])))
               names(regpara.summary)<-c(paste0('beta.',seq(1,n.beta) ),
                                      paste0('gam1.',seq(1,n.gamma1)),
                                      paste0('gam2.',seq(1,n.gamma2)))
               
              
       boot.reg<-data.frame(beta=t(regpara[1:n.beta]),gam1=t(regpara[mm1:mm2]),gam2=t(regpara[mm3:mm4]))
       
       boot.L1<-sapply(X=success.list,FUN=function(X) res1[[X]]$subdist.hazard1$Lambda)
       boot.L2<-sapply(X=success.list,FUN=function(X) res1[[X]]$subdist.hazard2$Lambda)
       boot.L1<-t(boot.L1)
       boot.L2<-t(boot.L2)
       time<-res1[[success.list[1]]]$subdist.hazard2$time
      
    }else{
      while(b.iter<n.boot){
        b.out<-bootstrap.result(Data)
        if(class(boot.out)=="list"){
          b.iter<-b.iter+1
          regpara<-as.numeric(as.character(b.out$reg.para$coef))
          regpara.summary<-data.frame(beta=t(regpara[1:n.beta]),gam1=t(regpara[mm1:mm2]),gam2=t(regpara[mm3:mm4]))
          L1<-b.out$subdist.hazard1$Lambda
          L2<-b.out$subdist.hazard2$Lambda
          row.names(L1)<-row.names(L2)<-NULL
          
          boot.reg<-rbind(boot.reg,regpara.summary)
          boot.L1<-rbind(boot.L1,L1)
          boot.L2<-rbind(boot.L2,L2)
          time<-b.out$subdist.hazard2$time
          remove(b.out)
      }#if
    }#while
    }#else
    
  output<-list()
  output$reg.coef<-boot.reg
  output$time<-time
  output$subdist.hazard1<-boot.L1
  output$subdist.hazard2<-boot.L2
  return(output)

}  # the end of the function ipw.bootstrap.ph2



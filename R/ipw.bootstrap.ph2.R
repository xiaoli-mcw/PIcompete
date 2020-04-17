#' To generate the empirical distribution for estimated sub-distribution hazard functions and regression coefficients, 
#' we adopted the two-phase weighted bootstrap procedure (Hyun et al, reference). We fit ipw.pi.competing function from main.R to 
#' the bootstrap weighted data.The most arguments except for three arguments, n.boot, parallel, n.core are necessary for fitting "ipw.pi.competing" function.
#' 
#' @importFrom parallel makeCluster clusterExport stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' 
#' @param n.boot the number of replications for bootstrap resampling 
#' @param parallel if parallel=TRUE, n.core should be specified. Default to FALSE.
#' @param n.core the nummber of cores for parallel computing
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
    boot.out<-ipw.pi.competing(Data=resamp.data,p.model=p.model,i.model1=i.model1,i.model2=i.model2,trans.r1=trans.r1,trans.r2=trans.r2,n.beta=n.beta, n.gamma1=n.gamma1,n.gamma2=n.gamma2,
                                   reg.initials=NULL,convergence.criteria=convergence.criteria,iteration.limit=iteration.limit,
                                   time.interval=time.interval,time.list=time.list,anal.var=FALSE)
                  
    return(boot.out)
  }
  ##################################################
  if(parallel&is.null(n.core)){
    stop('Specify the number of cores, "n.core" argument for parallel computing.')
  }
  
  
      
    if(parallel){
      cl <- makeCluster(n.core)
      clusterExport(cl,c("ipw.pi.competing"))
        
      registerDoParallel(cl)
      res1 <- foreach(n = 1:n.boot, .combine = c, .packages = c("fdrtool","nloptr"), .errorhandling=c('pass')) %dopar% 
        {b.out<-try(bootstrap.result(Data),silent = TRUE)
        b.out}
      
      reg.ind<-which(row.names(summary(res1))=="reg.coef")
      h1.ind<-which(row.names(summary(res1))=="subdist.hazard1")
      h2.ind<-which(row.names(summary(res1))=="subdist.hazard2")  
      
      regpara.summary<-sapply(X=reg.ind,FUN=function(X) as.numeric(as.character(res1[X]$reg.coef$coef)))
      boot.reg<-data.frame(t(regpara.summary))
      names(boot.reg)<-as.character(res1[reg.ind[1]]$reg.coef$label)
      
      L1<-sapply(X=h1.ind,FUN=function(X) res1[X]$subdist.hazard1$Lambda)
      time1<-res1[h1.ind[1]]$subdist.hazard1$time
      boot.L1<-t(L1)
        
      L2<-sapply(X=h2.ind,FUN=function(X) res1[X]$subdist.hazard2$Lambda)
      time2<-res1[h2.ind[1]]$subdist.hazard2$time
      boot.L2<-t(L2)
      if(length(reg.ind)<n.boot){
        n.fail<-n.boot-length(reg.ind)
        com<-paste0("Disconvergence of the algorithm arrises at ",n.fail," bootstrap samples")
       warning(com) 
      }
      stopCluster(cl) 
      
      } else{
        
        b.iter<-0
        boot.reg<-boot.L1<-boot.L2<-res1<-NULL
        
      while(b.iter<n.boot){
        b.out<-try(bootstrap.result(Data),silent = TRUE)
        if(class(b.out)=="list"){
          b.iter<-b.iter+1
          regpara<-as.numeric(as.character(b.out$reg.coef$coef))
          regpara.summary<-data.frame(t(regpara))
          
          L1<-b.out$subdist.hazard1$Lambda
          L2<-b.out$subdist.hazard2$Lambda
          row.names(L1)<-row.names(L2)<-NULL
          
          boot.reg<-rbind(boot.reg,regpara.summary)
          boot.L1<-rbind(boot.L1,L1)
          boot.L2<-rbind(boot.L2,L2)
          
          if(b.iter==n.boot){
            names(boot.reg)<-as.character(b.out$reg.coef$label)
            time1<-b.out$subdist.hazard1$time
            time2<-b.out$subdist.hazard2$time
          }
          remove(b.out)
      }#if
    }#while
    }#else
    
  output<-list()
  output$reg.coef<-boot.reg
  output$subdist.hazard1<-boot.L1
  output$subdist.hazard2<-boot.L2
  output$time1<-time1
  output$time2<-time2
  return(output)

}  # the end of the function ipw.bootstrap.ph2



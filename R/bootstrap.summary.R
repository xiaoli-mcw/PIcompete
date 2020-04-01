#' Based on the output obtained from ipw.bootstrap.ph2 function, empirical distribution based-variance/-confidence intervals for the estimated regression coefficients,
#' prevalences, cumulative sub-distribution hazards and cumulative incidences for events 1 and 2 are calculated. 
#' 
#' @import foreach
#'  
#' @param o.input the model fit output obtained from ipw.pi.competing function 
#' @param b.input the bootstrap output obtained from ipw.bootstrap.ph2 function
#' @param p.mat design matrix for predicting prevalence by using the inverse of logit function; 
#' both of vector and matrix types are allowed; the first component or column should include 1 for the intercept.
#' @param i.mat1 design matrix for predicting cumulative sub-distribution hazards and cumulative incidences for event 1
#' @param i.mat2 design matrix for predicting cumulative sub-distribution hazards and cumulative incidences for event 2
#' @param time.points time points at which cumulative sub-distribution hazards and cumulative incidences for events 1 and 2 are predicted
#' @param alpha The nominal coverage probability is (1-alpha)*100%. Default to 0.05
#' @return The output is a list of class ipw.bootstrap.ph2, which contains the following elements.
#'  \itemize{
#'  \item reg.covariance  empirical distribution based-covariance for the regression coefficients included in o.input
#'  \item reg.coef.ci empirical distribution based-confidence intervals for the regression coefficients included in o.input
#'  \item prev.ci empirical distribution based-confidence intervals for the predicted prevalence
#'  \item subdist.hazard1.ci confidence intervals for the predicted cumulative sub-distribution hazard function for event 1
#'  \item subdist.hazard2.ci confidence intervals for the predicted cumulative sub-distribution hazard function for event 2
#'  \item cum.inc1.ci confidence intervals for the predicted cumulative incidences for event 1
#'  \item cum.inc2.ci confidence intervals for the predicted cumulative incidences for event 2
#'  \item trans.r the transformation parameters used in the bootstrap summary
#'  \item p.mat the input information, the design matrix for prevalence
#'  \item i.mat1 the input information, the design matrix for cumulative sub-distribution hazard and cumulative incidence for event 1
#'  \item i.mat2 the input information, the design matrix for cumulative sub-distribution hazard and cumulative incidence for event 1
#'  

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

bootstrap.summary<-function(o.input,b.input,p.mat,i.mat1,i.mat2,time.points,alpha=0.05,...){

  ##########functions #######################################
  trans<-function(r,x){
    if(r==0){
      output<-x
    }else if(r>0){
      output<-log(1+r*x)/r
    }
    return(output)
  }
  
  smoothing<-function(time.list,L){
    s <- ksmooth(x=time.list, y=L, kernel = "box",bandwidth = 1,x.points=time.list)
    return(s$y)
  }
  
  risk.fn<-function(coef,mat){
    coef<-as.numeric(as.character(coef))
    if(class(mat)=="matrix"){
      risk<-mat%*%coef
    }else if(class(mat)=="numeric"){
      risk<-matrix(sum(mat*coef))
    }
    return(risk)
  }
  
  quantile.fn<-function(x,alpha){
    quantile(x,prob=c(alpha/2,1-alpha/2),na.rm=TRUE)
  }
  
  
  trans.hazard.fn<-function(time.list,subdist.hazard,risk,trans.r,time.points){
    L.e<-subdist.hazard*exp(risk)
    subdist.h<-trans(trans.r,L.e)
    cum.inc<-1-exp(-subdist.h)
    smooth.hazard<-smoothing(time.list,subdist.h)
    smooth.inc<-smoothing(time.list,cum.inc)
    
    hfun  <- stepfun(time.list,c(smooth.hazard[1],smooth.hazard), f = 1)
    ifun  <- stepfun(time.list,c(smooth.inc[1],smooth.inc), f = 1)
    #CADLAG function
    
    hazard<-hfun(time.points)
    incidence<-ifun(time.points)
    
    return(rbind(hazard,incidence))
  }
  ##### result from the original fitted model ###############
  
  ind1<-grep("logit",o.input$reg.coef$label)
  ind2<-grep("i.model1",o.input$reg.coef$label)
  ind3<-grep("i.model2",o.input$reg.coef$label)
  
  beta<-as.numeric(as.character(o.input$reg.coef$coef[ind1]))
  gam1<-as.numeric(as.character(o.input$reg.coef$coef[ind2]))
  gam2<-as.numeric(as.character(o.input$reg.coef$coef[ind3]))
  
  Lambda1<-o.input$subdist.hazard.fn1(time.points)
  Lambda2<-o.input$subdist.hazard.fn2(time.points)
  
  r1<-o.input$trans.r[1]
  r2<-o.input$trans.r[2]
  
  
  ############################
  
  if(class(p.mat)=="matrix"){
    n.design1<-nrow(p.mat)
  }else if(class(p.mat)=="numeric"){
    n.design1<-1
  }
  
  if(class(i.mat1)=="matrix"){
    n.design2<-nrow(i.mat1)
  }else if(class(i.mat1)=="numeric"){
    n.design2<-1
  }
  
  if(class(i.mat2)=="matrix"){
    n.design3<-nrow(i.mat2)
  }else if(class(i.mat2)=="numeric"){
    n.design3<-1
  }
  
  n.times<-length(time.points)
  
  ################ covariance matrix for regression coefficients estimates ###############
  b.means<-colMeans(b.input$reg.coef)
  reg.label<-colnames(b.input$reg.coef)
  n.regpara<-length(reg.label)
  
  b.rep<-nrow(b.input$reg.coef)
  b.mean.mat<-matrix(rep(b.means,b.rep),b.rep,n.regpara,2)
  
  c.reg.coef<-as.matrix(b.input$reg.coef-b.mean.mat)
  cov.mat<-(t(c.reg.coef)%*%c.reg.coef)/(b.rep-1)
  
  
  ############## empirical point-wise confidence intervals ###################
  ind1<-grep("logit",reg.label)
  ind2<-grep("i.model1",reg.label)
  ind3<-grep("i.model2",reg.label)
  
  
  p.risk<-sapply(X=1:b.rep,FUN=function(X) risk.fn(coef=b.input$reg.coef[X,ind1],mat=p.mat))
  i.risk1<-sapply(X=1:b.rep,FUN=function(X) risk.fn(coef=b.input$reg.coef[X,ind2],mat=i.mat1))
  i.risk2<-sapply(X=1:b.rep,FUN=function(X) risk.fn(coef=b.input$reg.coef[X,ind3],mat=i.mat2))

  op.risk<-risk.fn(coef=beta,mat=p.mat)
  oi.risk1<-risk.fn(coef=gam1,mat=i.mat1)  
  oi.risk2<-risk.fn(coef=gam2,mat=i.mat2)
  
  if(class(i.risk1)=="numeric"){
    i.risk1<-matrix(i.risk1,1,b.rep)
  }
  
  if(class(i.risk2)=="numeric"){
    i.risk2<-matrix(i.risk2,1,b.rep)
  }
  
  if(class(p.risk)=="numeric"){
    p.risk<-matrix(p.risk,1,b.rep)
  }
  colnames(i.risk2)<-colnames(i.risk1)<-colnames(p.risk)<-paste0(c("b.rep"),1:b.rep)
  
  ####### prevalence confidence interval ###################
  
  cp<-round((1-alpha)*100,1)
  label2<-paste0(cp,c("%_LL","%_UL"))
  
  ########## empirical confidence interval for regression coefficient estimates ########
   b.reg.mean<-colMeans(b.input$reg.coef)
   reg.mean<-b.reg.mean-as.numeric(as.character(o.input$reg.coef$coef))
   c.reg.coef<-b.input$reg.coef-matrix(rep(reg.mean,b.rep),b.rep,n.regpara,2)
   
   reg.coef.ci<-sapply(X=1:n.regpara,FUN=function(X) quantile.fn(c.reg.coef[,X],alpha))
   colnames(reg.coef.ci)<-names(b.input$reg.coef)
    
  ######### empirical confidence intervals for prevalence estimates ###############
  prev.bootstrap<-exp(p.risk)/(1+exp(p.risk))
  prev.mean<-rowMeans(prev.bootstrap)-exp(op.risk)/(1+exp(op.risk))
  c.prev.bootstrap<-prev.bootstrap-matrix(rep(prev.mean,each=b.rep),n.design1,b.rep,2)
  
  prev.ci<-sapply(X=1:n.design1,FUN=function(X) quantile.fn(c.prev.bootstrap[X,],alpha))
  rownames(prev.ci)<-label2
  colnames(prev.ci)<-paste0("p.mat row",seq(1,n.design1))
 
  ######### event 1: empirical confidence intervals for cumulative sub-distribution hazard and incidence estimates ###############
  o.L1.e<-foreach(x = 1:n.design2, .combine = rbind, .packages = c("foreach")) %do% {
    o.L1.e.out<-trans.hazard.fn(time.list=o.input$subdist.hazard1$time,subdist.hazard=o.input$subdist.hazard1$Lambda,
                                risk=oi.risk1[x,1],trans.r=r1,time.points)
    o.L1.e.out
  }   
  colnames(o.L1.e)<-time.points
  
  o.L1.e.h<-o.L1.e[which(rownames(o.L1.e)=="hazard"),]
  o.L1.e.i<-o.L1.e[which(rownames(o.L1.e)=="incidence"),]
  if(class(o.L1.e.h)=="numeric"){
    o.L1.e.h<-matrix(o.L1.e.h,1,n.times)
  }
  if(class(o.L1.e.i)=="numeric"){
    o.L1.e.i<-matrix(o.L1.e.i,1,n.times)
  }
  
  L1.e<-foreach(x = 1:n.design2, .combine = rbind, .packages = c("foreach")) %:%
         foreach(y=1:b.rep,.combine=rbind, .packages = c("foreach")) %do% {
    L1.e.out<-trans.hazard.fn(time.list=b.input$time1,subdist.hazard=b.input$subdist.hazard1[y,],
                     risk=i.risk1[x,y],trans.r=r1,time.points)
    L1.e.out
    } 
  colnames(L1.e)<-time.points
  h.ind<-which(rownames(L1.e)=="hazard" )
  i.ind<-which(rownames(L1.e)=="incidence")
  
  L1.e.h<-L1.e[h.ind,]
  L1.e.i<-L1.e[i.ind,]
  
  ############# to find quantiles ######################
  
  hazard1.q<-foreach(x = 1:n.design2, .combine = rbind, .packages = c("foreach")) %:%
    foreach(y=1:n.times,.combine=cbind, .packages = c("foreach")) %do% {
      shift1<-mean(L1.e.h[(1+(x-1)*b.rep):(x*b.rep),y])-o.L1.e.h[x,y]
      qout<- quantile.fn(L1.e.h[(1+(x-1)*b.rep):(x*b.rep),y]-shift1,alpha)
     qout
    }
  
  
  cuminc1.q<-foreach(x = 1:n.design2, .combine = rbind, .packages = c("foreach")) %:%
    foreach(y=1:n.times,.combine=cbind, .packages = c("foreach")) %do% {
      shift1<-mean(L1.e.i[(1+(x-1)*b.rep):(x*b.rep),y])-o.L1.e.i[x,y]
      qout<- quantile.fn(L1.e.i[(1+(x-1)*b.rep):(x*b.rep),y]-shift1,alpha)
      qout
    }
  
  ############# Labeling  ##########################
  hazard1.ci<-data.frame(cbind(Label=paste0("i.mat1 row",rep(seq(1,n.design2),each=2),  rep( paste0(", ",label2),n.design2)),hazard1.q))
  colnames(hazard1.ci)[-1]<-time.points
  rownames(hazard1.ci)<-NULL
   
  cuminc1.ci<-data.frame(cbind(Label=paste0("i.mat1 row",rep(seq(1,n.design2),each=2),  rep( paste0(", ",label2),n.design2)),cuminc1.q))
  colnames(cuminc1.ci)[-1]<-time.points
  rownames(cuminc1.ci)<-NULL
  
  
  
  ######### event 2: empirical confidence intervals for cumulative sub-distribution hazard and incidence estimates ###############
  o.L2.e<-foreach(x = 1:n.design3, .combine = rbind, .packages = c("foreach")) %do% {
      o.L2.e.out<-trans.hazard.fn(time.list=o.input$subdist.hazard2$time,subdist.hazard=o.input$subdist.hazard2$Lambda,
                                risk=oi.risk2[x,1],trans.r=r2,time.points)
      o.L2.e.out
    }   
  colnames(o.L2.e)<-time.points
  
  o.L2.e.h<-o.L2.e[which(rownames(o.L2.e)=="hazard"),]
  o.L2.e.i<-o.L2.e[which(rownames(o.L2.e)=="incidence"),]
  
  if(class(o.L2.e.h)=="numeric"){
    o.L2.e.h<-matrix(o.L2.e.h,1,n.times)
  }
  if(class(o.L2.e.i)=="numeric"){
    o.L2.e.i<-matrix(o.L2.e.i,1,n.times)
  }
  
  L2.e<-foreach(x = 1:n.design3, .combine = rbind, .packages = c("foreach")) %:%
    foreach(y=1:b.rep,.combine=rbind, .packages = c("foreach")) %do% {
      L2.e.out<-trans.hazard.fn(time.list=b.input$time1,subdist.hazard=b.input$subdist.hazard2[y,],
                                risk=i.risk2[x,y],trans.r=r2,time.points)
      L2.e.out
    } 
  colnames(L2.e)<-time.points
  h.ind2<-which(rownames(L2.e)=="hazard" )
  i.ind2<-which(rownames(L2.e)=="incidence")
  
  L2.e.h<-L2.e[h.ind2,]
  L2.e.i<-L2.e[i.ind2,]
  
  ############# to find quantiles ######################
  
  hazard2.q<-foreach(x = 1:n.design3, .combine = rbind, .packages = c("foreach")) %:%
    foreach(y=1:n.times,.combine=cbind, .packages = c("foreach")) %do% {
      shift2<-mean(L2.e.h[(1+(x-1)*b.rep):(x*b.rep),y])-o.L2.e.h[x,y]
      qout<- quantile.fn(L2.e.h[(1+(x-1)*b.rep):(x*b.rep),y]-shift2,alpha)
      qout
    }
  
  
  cuminc2.q<-foreach(x = 1:n.design3, .combine = rbind, .packages = c("foreach")) %:%
    foreach(y=1:n.times,.combine=cbind, .packages = c("foreach")) %do% {
      shift2<-mean(L2.e.i[(1+(x-1)*b.rep):(x*b.rep),y])-o.L2.e.i[x,y]
      qout<- quantile.fn(L2.e.i[(1+(x-1)*b.rep):(x*b.rep),y]-shift2,alpha)
      qout
    }
  
  ############# Labeling  ##########################
  hazard2.ci<-data.frame(cbind(Label=paste0("i.mat2 row",rep(seq(1,n.design3),each=2),  rep( paste0(", ",label2),n.design3)),hazard2.q))
  colnames(hazard2.ci)[-1]<-time.points
  rownames(hazard2.ci)<-NULL
  
  cuminc2.ci<-data.frame(cbind(Label=paste0("i.mat2 row",rep(seq(1,n.design3),each=2),  rep( paste0(", ",label2),n.design3)),cuminc2.q))
  colnames(cuminc2.ci)[-1]<-time.points
  rownames(cuminc2.ci)<-NULL
  
  
    
  output<-list()
  output$reg.covariance<-cov.mat
  output$reg.coef.ci<-reg.coef.ci
  output$prev.ci<-prev.ci
  output$subdist.hazard1.ci<-hazard1.ci
  output$subdist.hazard2.ci<-hazard2.ci
  output$cum.inc1.ci<-cuminc1.ci
  output$cum.inc2.ci<-cuminc2.ci
  output$trans.r<-c(trans.r1=r1,trans.r2=r2)
  output$p.mat<-p.mat
  output$i.mat1<-i.mat1
  output$i.mat2<-i.mat2
  
  
  return(output)
}
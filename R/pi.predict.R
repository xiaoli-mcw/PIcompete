#' Based on the output obtained by applying ipw.pi.competing function to data, 
#' pi.predict function calculates predicted prevalence by using the inverse function of logit function, 
#' cumulative sub-distribution hazards and cumulative incidences for events 1 and 2.
#' #' 
#' #@import 
#' 
#' @param input the output obtained from ipw.pi.competing 
#' @param p.mat design matrix for predicting prevalence by using the inverse of logit function; 
#' both of vector and matrix types are allowed; the first component or column should include 1 for the intercept.
#' @param i.mat1 design matrix for predicting cumulative sub-distribution hazards and cumulative incidences for event 1
#' @param i.mat2 design matrix for predicting cumulative sub-distribution hazards and cumulative incidences for event 2
#' @param time.points time points at which cumulative sub-distribution hazards and cumulative incidences for events 1 and 2 are predicted
#' @return The output is a list of class pi.predict, which contains the following elements.
#'  \itemize{
#'  \item prev predicted prevalence for the sub-groups with the covariates defined in p.mat
#'  \item subdist.hazard1 predicted cumulative sub-distribution hazards for event 1 for the sub-group with the covariates defined in i.mat1
#'  \item subdist.hazard2 predicted cumulative sub-distribution hazards for event 2 for the sub-group with the covariates defined in i.mat2
#'  \item cum.inc1 predicted cumulative incidences for event 1 for the sub-group with the covariates defined in i.mat1
#'  \item cum.inc2 predicted cumulative incidences for event 2 for the sub-group with the covariates defined in i.mat2
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

pi.pred<-function(input,p.mat,i.mat1,i.mat2,time.points){
  
  trans<-function(r,x){
    if(r==0){
      output<-x
    }else if(r>0){
      output<-log(1+r*x)/r
    }
    return(output)
  }
  
    ind1<-grep("logit",input$reg.coef$label)
    ind2<-grep("i.model1",input$reg.coef$label)
    ind3<-grep("i.model2",input$reg.coef$label)
    
    beta<-as.numeric(as.character(input$reg.coef$coef[ind1]))
    gam1<-as.numeric(as.character(input$reg.coef$coef[ind2]))
    gam2<-as.numeric(as.character(input$reg.coef$coef[ind3]))
    
    L1<-input$subdist.hazard.fn1(time.points)
    L2<-input$subdist.hazard.fn2(time.points)
    
    r1<-input$trans.r[1]
    r2<-input$trans.r[2]
  
  
  if(class(p.mat)=="matrix"){
    p.risk<-p.mat%*%beta
  }else if(class(p.mat)=="numeric"){
    p.risk<-sum(p.mat*beta)
  }
  
  if(class(i.mat1)=="matrix"){
    i.risk1<-i.mat1%*%gam1
  }else if(class(i.mat1)=="numeric"){
    i.risk1<-sum(i.mat1*gam1)
  }
  if(class(i.mat2)=="matrix"){
    i.risk2<-i.mat2%*%gam2
  }else if(class(i.mat2)=="numeric"){
    i.risk2<-sum(i.mat2*gam2)
  }
  
  
  if(length(time.points)>1){
    L1.e<-exp(i.risk1)%*%t(L1)
    L2.e<-exp(i.risk2)%*%t(L2)
  }else if(length(time.points)==1){
    L1.e<-exp(i.risk1)*L1
    L2.e<-exp(i.risk2)*L2
  }
  
  subdist.hazard1<-trans(r=r1,L1.e)
  subdist.hazard2<-trans(r=r2,L2.e)
  cum.inc1<-1-exp(-subdist.hazard1)
  cum.inc2<-1-exp(-subdist.hazard2)
  
  
  colnames(subdist.hazard1)<-colnames(subdist.hazard2)<-time.points
  colnames(cum.inc1)<-colnames(cum.inc2)<-time.points
  
  prev<-exp(p.risk)/(1+exp(p.risk))
  
  output<-list()
  output$prev<-prev
  output$subdist.hazard1<-subdist.hazard1
  output$subdist.hazard2<-subdist.hazard2
  output$cum.inc1<-cum.inc1
  output$cum.inc2<-cum.inc2
  return(output)
}


library(PIcompete)

trans.r1<-1
trans.r2<-1

model0<-"C ~ HR_status + risk_score"
model1<-"C + L1 + R1 ~ HR_status + risk_score"
model2<-"L2 + R2 ~ HR_status + risk_score"

n.beta<-3
n.gamma1<-2
n.gamma2<-2

true.reg.para<-c(-3.5,1,1,rep(0.3,2),-rep(0.3,2))

output <- ipw.pi.competing(Data=samp.data, p.model=model0, i.model1=model1, i.model2=model2, trans.r1=trans.r1, trans.r2=trans.r2, n.beta=n.beta, n.gamma1=n.gamma1, n.gamma2=n.gamma2)  #ipw.pi.competing

outputb <- ipw.bootstrap.ph2(n.boot=2, Data=samp.data, p.model=model0, i.model1=model1, i.model2=model2, trans.r1=trans.r1, trans.r2=trans.r2, n.beta=n.beta, n.gamma1=n.gamma1, n.gamma2=n.gamma2, iteration.limit=250, time.interval=0.1, time.list=seq(0,9,by=0.1))  #bootstrap

outputbp <- ipw.bootstrap.ph2(n.boot=2, Data=samp.data, p.model=model0, i.model1=model1, i.model2=model2, trans.r1=trans.r1, trans.r2=trans.r2, n.beta=n.beta, n.gamma1=n.gamma1, n.gamma2=n.gamma2, iteration.limit=250, time.interval=0.1, time.list=seq(0,9,by=0.1),parallel=TRUE,n.core=2)  #bootstrap parallel

x <- cbind(X0=1, samp.data[,c("HR_status","risk_score")])  #add 1 as the intercept col
x <- data.matrix(x)  #design matrix
time <- 10  #timepoints

sum1 <- bootstrap.summary(o.input = output, b.input = outputb, p.mat = x, i.mat1 = x[,-1], i.mat2 = x[,-1], time.points = time)  #summary
sum2 <- bootstrap.summary(o.input = output, b.input = outputbp, p.mat = x, i.mat1 = x[,-1], i.mat2 = x[,-1], time.points = time)  #summary


pred <- pi.pred(input = output, p.mat = x, i.mat1 = x[,-1], i.mat2 = x[,-1], time.points = time)  #predict

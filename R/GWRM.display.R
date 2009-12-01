GWRM.display<-function(model){
if (is.null(model$covars)==TRUE) model$covars<-"(Intercept)"
Table<-data.frame(covars=model$covars,estimates=model$betascoefs,se=model$se[1:length(model$betascoefs)],z=(model$betascoefs)/(model$se[1:length(model$betascoefs)]),p=2*pnorm(abs((model$betascoefs)/(model$se[1:length(model$betascoefs)])),lower=FALSE))
fit<-data.frame(loglikelihood=-model$optimum,AIC=model$aic,BIC=model$bic,df=model$df)
coefk<-model$coefficients[length(model$betascoefs)+1]
coefro<-model$coefficients[length(model$betascoefs)+2]
k<-model$betaIIpars[1]
ro<-model$betaIIpars[2]
betaII<-data.frame(par=c("k","ro"),coef=c(coefk,coefro),value=c(k,ro))
results<-list(Table=Table,betaII=betaII,Fit=fit,Convergence=model$code,Method=model$method)
return(results)
}
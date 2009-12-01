GWRM.fit<-function(formula,f=NULL,p0k=0,p0ro=1,p0beta=NULL,iters=10000,data,method=1){
#Design matrix
a <- model.frame(formula, data=data)
y<-model.extract(a, "response")
offset<-model.extract(a, "offset")
if (is.null(f)==TRUE) f<-rep(1,length(y))
if (is.null(offset)==TRUE){
	offset<-rep(1,length(y))
	covoffset<-FALSE}
	else covoffset<-TRUE
if (terms(formula)[[3]]==1){
	matrizmu<-matrix(1,c(length(y),1))
	namescovars<-c("(Intercept)")
		}
else {
  matrizmu<-model.matrix(terms(formula),model.frame(terms(formula),data=data,na.action=NULL))
  namescovars<-dimnames(matrizmu[0,])[[2]]
	}
ncovars<-ncol(matrizmu)
if (is.null(p0beta)==TRUE) p0beta<-rep(0,ncovars)

#Log-likelihood
logL<-function(p){
beta<-p[1:(ncovars)]
betak<-p[ncovars+1]
k<-exp(betak)
betaro<-p[ncovars+2]
ro<-1+exp(betaro)
mu<-offset*exp(matrizmu%*%beta)
a<-mu*(ro-1)/k
gama<-a+k+ro
-sum(f*(lgamma(a+y)-lgamma(a)+lgamma(k+y)-lgamma(k)-lgamma(gama+y)+lgamma(a+ro)+lgamma(k+ro)-lgamma(ro)))
}
#Optimizing log-likelihood
p0<-c(p0beta,p0k,p0ro)
if (method==1){
  fit<-nlm(logL,p=p0,hessian= TRUE,iterlim=iters)
  fit$value<-fit$minimum
  fit$par<-fit$estimate
  fit$convergence<-fit$code
  method="hessian nlm"
  }
if (method==2){
  fit<-optim(p0,logL,hessian=TRUE,control=list(maxit=iters))
  method<-"optim (hessian Nelder-Mead)"
  }
if (method==3){
  fit<-optim(p0,logL,method="BFGS",hessian=TRUE,control=list(maxit=iters))
  method<-"optim (BFGS)"
  }
if (method==4){
  fit<-optim(p0,logL,method="CG",hessian=TRUE,control=list(maxit=iters))
  method<-"optim (CG)"
  }
if (method==5){
  fit<-optim(p0,logL,method="L-BFGS-B",hessian=TRUE,control=list(maxit=iters))
  method<-"optim (L-BFGS-B)"
  }
if (method==6){
  fit<-optim(p0,logL,method="SANN",hessian=TRUE,control=list(maxit=iters))
  method<-"optim (SANN)"
  }
#Results
results<-list(
dataset=data.frame(data,offset=offset,f=f),
response=y,
model=formula,
covars=namescovars,
offset=covoffset,
optimum=fit$value+sum(f*lfactorial(y)),
aic=2*(fit$value+sum(f*lfactorial(y)))+(length(p0beta)+2)*2,
bic=2*(fit$value+sum(f*lfactorial(y)))+(length(p0beta)+2)*log(sum(f)),
df=sum(f)-(length(p0beta)+2),
coefficients=fit$par,
betaIIpars=c(exp(fit$par[ncovars+1]),1+exp(fit$par[ncovars+2])),
betascoefs=fit$par[1:ncovars],
fitted.values=offset*exp(matrizmu%*%fit$par[1:ncovars]),
hessian=fit$hessian,
cov=solve(fit$hessian),
se=sqrt(diag(solve(fit$hessian))),
corr=solve(fit$hessian)/(sqrt(diag(solve(fit$hessian)))%o%sqrt(diag(solve(fit$hessian)))),
code=fit$convergence,
method=method
)
return(results)
}
GWRM.stats<-function(model=NULL,covars=NULL,alpha=0.05){
if (is.null(covars)==TRUE){ 
	niveles<-1
	covars<-"No"
	matrizmu<-1
  offset<-1}
else {
  ncasos<-nrow(covars)
  niveles<-data.frame(rep(1,nrow(covars)),covars)
  names(niveles)[1]<-as.character(model$model[[2]])
  matrizmu<-model.matrix(terms(model$model),model.frame(terms(model$model),data=niveles,na.action=NULL))
  a <- model.frame(model$model, data=niveles)
  offset<-model.extract(a, "offset")
  }
if (is.null(offset)==TRUE){
	offset<-rep(1,nrow(covars))}  
ncovars<-ncol(matrizmu)
beta<-model$coefficients[1:(ncovars)]
k<-exp(model$coefficients[ncovars+1])
ro<-1+exp(model$coefficients[ncovars+2])
mus<-offset*exp(matrizmu%*%beta)
liminf<-vector(mode="numeric",length=ncasos)
limsup<-vector(mode="numeric",length=ncasos)
for (i in 1:ncasos)
{error<-mus[i]*sqrt(matrizmu[i,]%*%model$cov[1:ncovars,1:ncovars]%*%matrizmu[i,])
liminf[i]<-mus[i]-qnorm(1-alpha/2)*error
limsup[i]<-mus[i]+qnorm(1-alpha/2)*error}
a<-mus*(ro-1)/k
var<-mus*((a+ro-1)*(k+ro-1))/((ro-1)*(ro-2))
prand<-mus/var
pliabi<-((ro-1)*(k+1))/((a+ro-1)*(k+ro-1))
pprone<-a/(a+ro-1)
results<-list(
params=data.frame(Levels=covars,a=a,k=k,ro=ro),
stats=data.frame(Levels=covars,Mean.est=mus,Low.bound=liminf,Upp.bound=limsup,Var.est=var),partvar=data.frame(Levels=covars,Randomness=prand,Liability=pliabi,Proneness=pprone))
return(results)
}
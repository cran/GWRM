GWRM.add<-function(object=NULL,scope=NULL,method=1,iters=10000){
dataset<-object$dataset
vblesfuera<-terms.formula(scope)
vblesdentro<-terms.formula(object$model)
niters<-length(attributes(vblesfuera)$term.labels)
covs<-c();chis<-c();dfs<-c();ps<-c()
for (k in 1:niters){
	entra<-vblesfuera[k]
  modelo1<-GWRM.fit(formula=update.formula(object$model,as.formula(paste("~ . +", terms(entra)[[2]])),evaluate = FALSE),f=dataset$f,data=dataset,iters=iters,method=method)
  covs[k]<-attributes(vblesfuera)$term.labels[k]
	chis[k]<-2*(object$optimum-modelo1$optimum)
	dfs[k]<-length(modelo1$coefficients)-length(object$coefficients)
	ps[k]<-pchisq(2*(object$optimum-modelo1$optimum),length(modelo1$coefficients)-length(object$coefficients),lower.tail=FALSE)
	}
data.frame(Covariates=covs,chisq=chis,df=dfs,p=ps)
}
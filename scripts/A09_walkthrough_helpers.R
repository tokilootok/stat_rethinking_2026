# Shared numerical helpers for A09 problem walkthroughs. Base R only.
# Grid quadrature approximates a posterior; it is not MCMC or an exact symbolic fit.
a09_normalize <- function(z) {w<-exp(z-max(z));w/sum(w)}
a09_logsum <- function(z) {m<-max(z);m+log(sum(exp(z-m)))}
a09_summary <- function(x) c(mean=mean(x),lo90=unname(quantile(x,.05)),
  hi90=unname(quantile(x,.95)))
a09_draw <- function(f,S=20000) {
  i<-sample.int(nrow(f$parameters),S,replace=TRUE,prob=f$mass)
  f$parameters[i,,drop=FALSE]
}
a09_regression <- function(x,y,family="binomial",size=rep(1,length(y)),
  offset=rep(0,length(y)),prior_mean=c(0,0),prior_sd=c(1,1),
  agrid=seq(-7,7,length.out=201),bgrid=seq(-6,6,length.out=201)) {
  stopifnot(length(x)==length(y),length(size)==length(y),length(offset)==length(y))
  pars<-expand.grid(alpha=agrid,beta=bgrid)
  logprior<-dnorm(pars$alpha,prior_mean[1],prior_sd[1],log=TRUE)+
    dnorm(pars$beta,prior_mean[2],prior_sd[2],log=TRUE)
  ll<-matrix(0,nrow(pars),length(y))
  for(i in seq_along(y)) {
    eta<-pars$alpha+pars$beta*x[i]+offset[i]
    ll[,i]<-if(family=="binomial") dbinom(y[i],size[i],plogis(eta),log=TRUE)
      else dpois(y[i],exp(eta),log=TRUE)
  }
  logmass<-logprior+rowSums(ll);mass<-a09_normalize(logmass)
  edge<-pars$alpha %in% range(agrid)
  if(length(bgrid)>1) edge<-edge | pars$beta %in% range(bgrid)
  stopifnot(sum(mass[edge])<1e-5)
  list(parameters=pars,mass=mass,logmass=logmass,loglik=ll,edge=edge)
}
a09_cells <- function(y,n,family="binomial",prior_mean=0,prior_sd=1,
                      S=20000,step=.01) {
  a<-seq(-12,8,by=step)
  likelihood_mean<-if(family=="binomial") plogis(a) else exp(a)
  draws<-matrix(NA_real_,S,length(y));means<-numeric(length(y))
  for(k in seq_along(y)) {
    ll<-if(family=="binomial") dbinom(y[k],n[k],likelihood_mean,log=TRUE)
      else dpois(y[k],n[k]*likelihood_mean,log=TRUE)
    w<-a09_normalize(ll+dnorm(a,prior_mean,prior_sd,log=TRUE))
    stopifnot(sum(w[a < -11 | a > 7])<1e-6)
    means[k]<-sum(w*likelihood_mean)
    draws[,k]<-sample(likelihood_mean,S,replace=TRUE,prob=w)
  }
  list(draws=draws,mean=means)
}
a09_loo <- function(f) {
  # Recompute the grid posterior without each observation, then integrate its likelihood.
  # This is numerical leave-one-out integration, not PSIS on posterior samples.
  sapply(seq_len(ncol(f$loglik)),function(i) {
    train<-f$logmass-f$loglik[,i]
    stopifnot(sum(a09_normalize(train)[f$edge])<1e-4)
    a09_logsum(train+f$loglik[,i])-a09_logsum(train)
  })
}
a09_output <- function(problem) {
  path<-file.path("courses","figures",paste0("A09_problem",problem))
  dir.create(path,recursive=TRUE,showWarnings=FALSE);path
}

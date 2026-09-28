# Beginner demonstrations for A09 sections 39, 65, 74, 77.
# Run from repository root: Rscript scripts/A09_lab_walkthroughs.R
# Requires posterior for the diagnostic demonstration; no grant model is fitted.

cat("BEGIN_counts\n")
mu <- 4; phi <- 2; zi <- .30; lambda <- 4
moments <- data.frame(model=c("Poisson", "Gamma-Poisson", "Zero-inflated Poisson"),
  mean=c(mu,mu,(1-zi)*lambda),
  variance=c(mu,mu+mu^2/phi,(1-zi)*lambda+zi*(1-zi)*lambda^2),
  prob_zero=c(dpois(0,mu),dnbinom(0,size=phi,mu=mu),zi+(1-zi)*dpois(0,lambda)))
print(moments,digits=4,row.names=FALSE)
# Independent latent intensities yield extra-Poisson variation.
set.seed(9039)
latent_rate <- rgamma(100000,shape=phi,rate=phi/mu)
y <- rpois(length(latent_rate),latent_rate)
cat("Simulated Gamma-Poisson mean and variance:",round(mean(y),3),round(var(y),3),"\n")
stopifnot(abs(mean(y)-mu)<.1,abs(var(y)-(mu+mu^2/phi))<.5)
# The installed rethinking helper uses scale = mu / phi.
if (requireNamespace("rethinking",quietly=TRUE))
  stopifnot(max(abs(rethinking::dgampois(0:30,mu,scale=mu/phi)-
                    dnbinom(0:30,size=phi,mu=mu)))<1e-12)
cat("END_counts\n")

cat("BEGIN_diagnostics\n")
# Synthetic chains, not the NWOGrants fit; posterior supplies modern diagnostics.
set.seed(9065)
chains <- matrix(rnorm(4000),nrow=1000,ncol=4)
shifted <- chains
shifted[,4] <- shifted[,4]+2
slow <- replicate(4,as.numeric(arima.sim(list(ar=.95),n=1000,
                                       sd=sqrt(1-.95^2))))
report <- function(x) c(Rhat=posterior::rhat(x),
  ESS_bulk=posterior::ess_bulk(x),ESS_tail=posterior::ess_tail(x),
  MCSE_mean=posterior::mcse_mean(x))
print(round(rbind(Independent=report(chains),Shifted_chain=report(shifted),
                  Correlated=report(slow)),3))
cat("If posterior SD = 0.03 and ESS for the mean = 100, MCSE is about",
    .03/sqrt(100),"\n")
cat("At ESS = 2500, it is about",.03/sqrt(2500),"\n")
cat("END_diagnostics\n")

cat("BEGIN_sbc\n")
# Exact Beta-binomial toy model: distinct from the grants' Normal-logit prior.
# One fixed numerical illustration of the rank calculation:
p_true <- .30; illustrative_draws <- c(.10,.20,.35,.45)
cat("Illustrative rank:",sum(illustrative_draws<p_true),"out of 0:4\n")
# If 6 of 20 trials succeed, Beta(2,2) updates to Beta(8,16).
cat("Posterior mean for y=6, n=20:",8/(8+16),"\n")
set.seed(9074)
B <- 1000; L <- 19; n <- 20
ranks <- biased_ranks <- integer(B)
for (b in seq_len(B)) {
  truth <- rbeta(1,2,2)                    # 1: draw a fresh truth from the prior
  y <- rbinom(1,n,truth)                  # 2: simulate its data
  draws <- rbeta(L,2+y,2+n-y)             # 3: exact independent posterior draws
  ranks[b] <- sum(draws<truth)            # 4: number below truth, from 0 to 19
  wrong_draws <- plogis(qlogis(draws)+.8) # deliberately biased algorithm
  biased_ranks[b] <- sum(wrong_draws<truth)
}
bins <- function(r) tabulate(r %/% 4+1,nbins=5)
print(data.frame(rank_range=c("0-3","4-7","8-11","12-15","16-19"),
  expected=B/5,exact=bins(ranks),biased=bins(biased_ranks)),row.names=FALSE)
dir.create("courses/figures/A09_lab",recursive=TRUE,showWarnings=FALSE)
png("courses/figures/A09_lab/sbc.png",width=1200,height=600,res=130)
par(mfrow=c(1,2),mar=c(4,4,3,1))
for (r in list(ranks,biased_ranks)) {
  barplot(bins(r),names.arg=c("0-3","4-7","8-11","12-15","16-19"),
    ylim=c(0,max(c(bins(ranks),bins(biased_ranks)))+20),
    main=if(identical(r,ranks)) "Exact posterior" else "Deliberately shifted posterior",
    xlab="Rank among 19 draws (grouped)",ylab="Number of simulations",col="steelblue")
  abline(h=B/5,lty=2)
}
invisible(dev.off())
cat("END_sbc\n")

cat("BEGIN_pooling\n")
# Teaching approximation with known hyperparameters, not a fitted grant model.
m <- 0; tau <- .30
b_hat <- c(-.80,-.80); se <- c(.10,.80)
w <- tau^2/(tau^2+se^2)
pooled <- w*b_hat+(1-w)*m
print(data.frame(discipline=c("Precise","Sparse"),estimate=b_hat,se=se,
  data_weight=round(w,3),pooled_mean=round(pooled,3)),row.names=FALSE)
# A bivariate Normal prior for baseline and contrast:
tau_alpha <- .50; tau_beta <- .30; rho <- .40
Sigma <- matrix(c(tau_alpha^2,rho*tau_alpha*tau_beta,
                  rho*tau_alpha*tau_beta,tau_beta^2),2,2)
print(Sigma)
set.seed(9077)
z <- matrix(rnorm(20000),ncol=2)
coef <- sweep(z %*% chol(Sigma),2,c(qlogis(.20),0),FUN="+")
cat("Simulated prior correlation:",round(cor(coef)[1,2],3),"\n")
stopifnot(max(abs(cov(coef)-Sigma))<.02)
png("courses/figures/A09_lab/pooling.png",width=1000,height=600,res=130)
par(mar=c(5,5,3,1))
plot(b_hat,1:2,xlim=c(-.9,.1),ylim=c(.5,2.5),yaxt="n",pch=1,
  xlab="Discipline log-odds contrast",ylab="",main="Partial pooling: a normal approximation")
axis(2,at=1:2,labels=c("Precise","Sparse"),las=1)
arrows(b_hat,1:2,pooled,1:2,length=.1,col="steelblue",lwd=2)
points(pooled,1:2,pch=19,col="steelblue"); abline(v=m,lty=2)
legend("bottomright",c("Unpooled estimate","Pooled mean","Population mean"),
  pch=c(1,19,NA),lty=c(NA,NA,2),col=c("black","steelblue","black"),bty="n")
invisible(dev.off())
cat("END_pooling\n")

# Reproduce the tables and figures in A09 sections 8 and 9.
# Run from the repository root: Rscript --vanilla scripts/A09_prior_moments.R
# Base R only; no posterior fitting or external data are required.
out_dir <- "courses/figures/A09_priors"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Normal(mu, sigma) uses standard deviation sigma.
# Integrating over Z ~ Normal(0,1) avoids endpoint singularities on p in (0,1).
logistic_normal_moments <- function(mu, sigma) {
  expect <- function(f) integrate(function(z) f(plogis(mu + sigma*z))*dnorm(z),
                                  -Inf, Inf, rel.tol=1e-10, subdivisions=1000L)$value
  m <- expect(identity)
  v <- expect(function(p) (p-m)^2)
  c(mean=m, median=plogis(mu), variance=v, sd=sqrt(v),
    skewness=expect(function(p) (p-m)^3)/v^(3/2),
    excess_kurtosis=expect(function(p) (p-m)^4)/v^2-3,
    q025=plogis(mu+qnorm(.025)*sigma), q975=plogis(mu+qnorm(.975)*sigma),
    extreme=pnorm(qlogis(.01),mu,sigma)+pnorm(qlogis(.99),mu,sigma,lower.tail=FALSE))
}
beta_moments <- function(a,b) {
  v <- a*b/((a+b)^2*(a+b+1))
  c(mean=a/(a+b), median=qbeta(.5,a,b), variance=v, sd=sqrt(v),
    skewness=2*(b-a)*sqrt(a+b+1)/((a+b+2)*sqrt(a*b)),
    excess_kurtosis=6*((a-b)^2*(a+b+1)-a*b*(a+b+2))/(a*b*(a+b+2)*(a+b+3)),
    q025=qbeta(.025,a,b),q975=qbeta(.975,a,b),
    extreme=pbeta(.01,a,b)+pbeta(.99,a,b,lower.tail=FALSE))
}
intercepts <- rbind(
  logistic_normal_moments(0,.5), logistic_normal_moments(0,1),
  logistic_normal_moments(0,2.5), logistic_normal_moments(0,10),
  logistic_normal_moments(qlogis(.1),1),
  beta_moments(1,1), beta_moments(2,2), beta_moments(2,18))
intercepts <- data.frame(prior=c("alpha ~ Normal(0, 0.5)","alpha ~ Normal(0, 1)",
  "alpha ~ Normal(0, 2.5)","alpha ~ Normal(0, 10)",
  "alpha ~ Normal(logit(0.1), 1)","p0 ~ Beta(1, 1)","p0 ~ Beta(2, 2)","p0 ~ Beta(2, 18)"),intercepts)
intercepts$count_variance_N20 <- with(intercepts,20*mean*(1-mean)+20*19*variance)
write.csv(intercepts,file.path(out_dir,"intercept_moments.csv"),row.names=FALSE)

slope_sd <- c(.25,.5,1,2)
odds <- data.frame(sigma_beta=slope_sd, mean_beta=0, variance_beta=slope_sd^2,
  median_OR=1,mean_OR=exp(slope_sd^2/2),
  variance_OR=expm1(slope_sd^2)*exp(slope_sd^2),
  q025_OR=exp(qnorm(.025)*slope_sd),q975_OR=exp(qnorm(.975)*slope_sd))
write.csv(odds,file.path(out_dir,"slope_odds_moments.csv"),row.names=FALSE)
# Independent alpha ~ Normal(logit(.1), .5), beta ~ Normal(0, sigma_beta).
# For fixed x, eta is normal, so the same one-dimensional integration applies.
probability_rows <- do.call(rbind,lapply(slope_sd,function(sb) {
  do.call(rbind,lapply(c(0,1,2),function(x) {
    data.frame(sigma_beta=sb,x=x,
      as.list(logistic_normal_moments(qlogis(.1),sqrt(.5^2+x^2*sb^2))))
  }))
}))
write.csv(probability_rows,file.path(out_dir,"slope_probability_moments.csv"),row.names=FALSE)

# Figure 1: same prior represented before and after the inverse link.
cols <- c("#0072B2","#009E73","#D55E00","#CC79A7")
sds <- c(.5,1,2.5,10)
png(file.path(out_dir,"intercept_scales.png"),width=1800,height=800,res=150)
par(mfrow=c(1,2),mar=c(4.5,4.5,3,1),las=1)
a <- seq(-12,12,length.out=1201)
matplot(a,sapply(sds,function(z) dnorm(a,0,z)),type="l",lty=1,lwd=2,
        col=cols,xlab="Intercept alpha (log odds)",ylab="Density",
        main="Normal priors on the intercept")
legend("topright",legend=paste("SD =",sds),col=cols,lty=1,lwd=2,bty="n",cex=.85)
breaks <- seq(0,1,by=.02)
bin_mass <- sapply(sds,function(z) diff(pnorm(qlogis(breaks),0,z)))
matplot(head(breaks,-1)+.01,bin_mass,type="l",lty=1,lwd=2,col=cols,
        xlab="Baseline probability p0",ylab="Prior probability per 0.02-wide bin",
        main="Same priors after inverse logit",xlim=c(0,1))
legend("top",legend=paste("SD =",sds),col=cols,lty=1,lwd=2,bty="n",cex=.85)
dev.off()

# Figure 2: moments as the intercept SD changes; same E[p0] but different spread.
ss <- seq(.05,10,length.out=140)
v <- vapply(ss,function(z) logistic_normal_moments(0,z)["variance"],numeric(1))
png(file.path(out_dir,"intercept_variance.png"),width=1800,height=750,res=150)
par(mfrow=c(1,2),mar=c(4.5,4.5,3,1),las=1)
plot(ss,v,type="l",lwd=3,col=cols[1],ylim=c(0,.25),
     xlab="Intercept prior SD",ylab="Variance of p0",main="Probability variance approaches 0.25")
abline(h=1/12,col="gray45",lty=2)
legend("bottomright",c("Logistic-normal","Uniform(0,1): variance 1/12"),
       col=c(cols[1],"gray45"),lty=c(1,2),lwd=c(3,1),bty="n",cex=.8)
plot(ss,5+380*v,type="l",lwd=3,col=cols[3],ylim=c(0,100),
     xlab="Intercept prior SD",ylab="Prior predictive variance of K",
     main="K | p0 ~ Binomial(20, p0)")
abline(h=5,col="gray45",lty=2)
legend("bottomright",c("Uncertain shared p0","Fixed p0 = 0.5: variance 5"),
       col=c(cols[3],"gray45"),lty=c(1,2),lwd=c(3,1),bty="n",cex=.8)
dev.off()

# Figure 3: share random standard-normal draws across panels for comparison.
set.seed(20260923)
xgrid <- seq(-2,2,length.out=121)
a_draw <- rnorm(40,qlogis(.1),.5)
z_beta <- rnorm(40)
png(file.path(out_dir,"slope_probability_curves.png"),width=1700,height=1250,res=150)
par(mfrow=c(2,2),mar=c(4,4,3,1),las=1)
for (k in seq_along(slope_sd)) {
  pp <- plogis(outer(xgrid,z_beta*slope_sd[k])+matrix(a_draw,length(xgrid),40,byrow=TRUE))
  matplot(xgrid,pp,type="l",lty=1,col=adjustcolor(cols[k],alpha.f=.3),ylim=c(0,1),
          xlab="Standardized predictor x",ylab="Probability p(x)",
          main=paste("beta ~ Normal(0,",slope_sd[k],")"))
  abline(v=0,col="gray70",lty=3)
  lines(xgrid,rep(.1,length(xgrid)),lty=2,lwd=2)
}
dev.off()

# Numerical checks: symmetry, boundedness, analytic beta benchmarks and x=0 invariance.
stopifnot(max(abs(intercepts$mean[1:4]-.5))<1e-8,
          all(intercepts$variance>0 & intercepts$variance<=.25),
          abs(intercepts$variance[6]-1/12)<1e-12,
          abs(intercepts$excess_kurtosis[6]+1.2)<1e-12)
at_zero <- probability_rows[probability_rows$x==0,]
stopifnot(diff(range(at_zero$mean))<1e-12,diff(range(at_zero$variance))<1e-12)
# Independent simulation check of the asymmetric case used in the text.
set.seed(17)
p_check <- plogis(rnorm(300000,qlogis(.1),1))
stopifnot(abs(mean(p_check)-intercepts$mean[5])<.002,
          abs(var(p_check)-intercepts$variance[5])<.001)
cat("Wrote 3 figures and 3 moment tables to",out_dir,"\n")
print(intercepts,digits=4,row.names=FALSE)
print(odds,digits=4,row.names=FALSE)

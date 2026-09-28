# A10 worked examples: run Rscript scripts/A10_walkthroughs.R from repository root.
# Base R, except posterior for the synthetic MCMC diagnostic demonstration.
# Fits fixed Gaussian-latent sensitivity scenarios by deterministic quadrature.

# Section 6
cat("BEGIN_6\n")
# Same data-generating probabilities as scripts/A10_sensitivity.R.
# Base-R Bernoulli draws; this is simulated data, not UCBadmit.
set.seed(12)
N <- 2000
G <- sample(1:2,N,replace=TRUE)
U <- rbinom(N,1,.1)
D <- 1+rbinom(N,1,ifelse(G==1,U,.75))
p0 <- matrix(c(.1,.1,.1,.3),nrow=2) # rows: D; columns: G
p1 <- matrix(c(.3,.3,.5,.5),nrow=2)
p <- ifelse(U==1,p1[cbind(D,G)],p0[cbind(D,G)])
A <- rbinom(N,1,p)
sim <- data.frame(G,D,U,A)
counts <- aggregate(cbind(admitted=A,applied=rep(1,N))~G+D,sim,sum)
counts$rate <- with(counts,admitted/applied)
print(counts,row.names=FALSE)
# Population truths, obtained by averaging over Pr(U=1)=0.1.
truth <- data.frame(department=1:2,
  U_high_G1=c(0,1),U_high_G2=c(.1,.1),
  observed_G1=c(.1,.3),observed_G2=c(.14,.32),
  intervention_G1=c(.12,.12),intervention_G2=c(.14,.32))
truth$observed_difference <- with(truth,observed_G1-observed_G2)
truth$controlled_difference <- with(truth,intervention_G1-intervention_G2)
print(truth,row.names=FALSE)
cat("True total contrast G1-G2:",.12-(.25*.14+.75*.32),"\n")
stopifnot(all.equal(truth$controlled_difference,c(-.02,-.20)))
cat("END_6\n")

# Section 8
cat("BEGIN_8\n")
summary90 <- function(x) c(mean=mean(x),lo90=unname(quantile(x,.05)),
                           hi90=unname(quantile(x,.95)))
normalize <- function(z) {w<-exp(z-max(z)); w/sum(w)}
set.seed(1008)
S <- 30000
agrid <- seq(-8,8,length.out=401)
totals <- aggregate(cbind(admitted,applied)~G,counts,sum)
naive_p <- sapply(1:2,function(g) {
  w <- normalize(dbinom(totals$admitted[g],totals$applied[g],
                       plogis(agrid),log=TRUE)+dnorm(agrid,log=TRUE))
  sample(plogis(agrid),S,replace=TRUE,prob=w)
})
cat("Naive aggregate posterior, G1-G2, percentage points:\n")
print(round(100*summary90(naive_p[,1]-naive_p[,2]),2))
cat("END_8\n")

# Section 9
cat("BEGIN_9\n")
set.seed(1009)
naive_cell <- sapply(seq_len(nrow(counts)),function(k) {
  w <- normalize(dbinom(counts$admitted[k],counts$applied[k],
                       plogis(agrid),log=TRUE)+dnorm(agrid,log=TRUE))
  sample(plogis(agrid),S,replace=TRUE,prob=w)
})
cell_index <- function(g,d) which(counts$G==g & counts$D==d)
naive_difference <- sapply(1:2,function(d)
  naive_cell[,cell_index(1,d)]-naive_cell[,cell_index(2,d)])
colnames(naive_difference) <- c("Department_1","Department_2")
print(round(100*t(apply(naive_difference,2,summary90)),2))
cat("END_9\n")

# Section 12
cat("BEGIN_12\n")
u <- c(-1,0,1)
# Illustrative baselines, not fitted coefficients: admission .1, selection .25.
print(data.frame(u=u,
  admission=round(plogis(qlogis(.1)+1*u),3),
  dept2_group1=round(plogis(qlogis(.25)+1*u),3),
  dept2_group2=round(plogis(qlogis(.25)+0*u),3)),row.names=FALSE)
cat("Odds multiplier for coefficient 1:",round(exp(1),3),"\n")
cat("END_12\n")

# Section 14
cat("BEGIN_14\n")
# Uniform log-odds coefficient does not imply uniform odds ratio.
z <- c(.05,.50,.95)
print(data.frame(prior_quantile=z,coefficient=z,odds_ratio=round(exp(z),3),
  probability_from_10pct=round(plogis(qlogis(.1)+z),3)),row.names=FALSE)
cat("Prior mean odds ratio:",round(exp(1)-1,3),"\n")
cat("END_14\n")

# Section 17
cat("BEGIN_17\n")
# Fixed b,g scenarios: integrate out each new applicant's u ~ Normal(0,1).
# Conditional on delta, the two department intercepts are independent.
fit_group <- function(group,b,g,K=201,Q=301,S=30000) {
  aa <- seq(-8,8,length.out=K)
  dd <- aa
  uu <- qnorm(((1:Q)-.5)/Q) # equal-weight quadrature under standard Normal
  admission <- plogis(outer(aa,b*uu,"+"))
  selection <- plogis(outer(uu*g,dd,"+"))
  qbar <- colMeans(selection)
  # P(A=1 | G,D,parameters), with selection-dependent latent weights.
  joint2 <- admission %*% selection/Q
  joint1 <- admission %*% (1-selection)/Q
  p2 <- sweep(joint2,2,qbar,"/")
  p1 <- sweep(joint1,2,1-qbar,"/")
  cells <- counts[counts$G==group,]
  cells <- cells[match(1:2,cells$D),]
  lw <- lapply(list(p1,p2),function(p) matrix(0,nrow=K,ncol=K))
  probs <- list(p1,p2)
  for (d in 1:2) lw[[d]] <- matrix(dbinom(cells$admitted[d],
    cells$applied[d],as.vector(probs[[d]]),log=TRUE),K,K)+dnorm(aa,log=TRUE)
  logsum <- function(x) {v<-max(x);v+log(sum(exp(x-v)))}
  z <- lapply(lw,function(x) apply(x,2,logsum))
  delta_mass <- normalize(dbinom(cells$applied[2],sum(cells$applied),
    qbar,log=TRUE)+dnorm(dd,log=TRUE)+z[[1]]+z[[2]])
  amass <- lapply(lw,function(x) apply(x,2,normalize))
  stopifnot(sum(delta_mass[abs(dd)>7])<1e-6)
  for (d in 1:2) stopifnot(sum(colSums(amass[[d]][abs(aa)>7,,drop=FALSE])*
                             delta_mass)<1e-6)
  # Population-standardized admission probability at each intercept.
  target_p <- rowMeans(admission)
  exact_mean <- sapply(1:2,function(d)
    sum(delta_mass*colSums(amass[[d]]*target_p)))
  di <- sample.int(K,S,replace=TRUE,prob=delta_mass)
  ai <- matrix(NA_integer_,S,2)
  for (j in unique(di)) for (d in 1:2) {
    rows <- which(di==j)
    ai[rows,d] <- sample.int(K,length(rows),replace=TRUE,prob=amass[[d]][,j])
  }
  list(mean=exact_mean,a=matrix(aa[ai],S,2),delta=dd[di],b=b,g=g,
       standardized=matrix(target_p[ai],S,2))
}
cat("Model: independent Normal(0,1) priors on a and delta; fixed b,g.\n")
cat("END_17\n")

# Section 19
cat("BEGIN_19\n")
settings <- list(No_latent=c(0,0,0,0),
  Outcome_only=c(1,1,0,0),Selection_only=c(0,0,1,0),
  Official_fixed=c(1,1,1,0),Stronger=c(2,2,2,0)) # b1,b2,g1,g2
set.seed(1019)
fits <- lapply(settings,function(v) list(
  fit_group(1,v[1],v[3]),fit_group(2,v[2],v[4])))
scenario_draws <- lapply(fits,function(f) f[[1]]$standardized-f[[2]]$standardized)
scenario_table <- do.call(rbind,lapply(names(settings),function(nm) {
  x<-scenario_draws[[nm]]
  data.frame(scenario=nm,department=1:2,
    round(100*t(apply(x,2,summary90)),2),prob_positive=round(colMeans(x>0),3))
}))
print(scenario_table,row.names=FALSE)
# Deterministic posterior means, not random draws, for resolution check.
v <- settings$Official_fixed
fine <- list(fit_group(1,v[1],v[3],K=401,Q=601,S=10),
             fit_group(2,v[2],v[4],K=401,Q=601,S=10))
coarse_mean <- fits$Official_fixed[[1]]$mean-fits$Official_fixed[[2]]$mean
fine_mean <- fine[[1]]$mean-fine[[2]]$mean
cat("Maximum refinement change in contrast (percentage points):",
    round(100*max(abs(coarse_mean-fine_mean)),4),"\n")
stopifnot(max(abs(coarse_mean-fine_mean))<.001)
# Null latent case must agree with independent cell-logit grid integration.
null_mean <- fits$No_latent[[1]]$mean-fits$No_latent[[2]]$mean
independent_mean <- sapply(1:2,function(d) {
  m <- sapply(1:2,function(g) {
    k<-cell_index(g,d)
    w<-normalize(dbinom(counts$admitted[k],counts$applied[k],plogis(agrid),
                       log=TRUE)+dnorm(agrid,log=TRUE))
    sum(w*plogis(agrid))
  });m[1]-m[2]
})
stopifnot(max(abs(null_mean-independent_mean))<1e-5)
dir.create("courses/figures/A10",recursive=TRUE,showWarnings=FALSE)
png("courses/figures/A10/scenarios.png",width=1300,height=650,res=130)
par(mfrow=c(1,2),mar=c(5,9,3,1))
for(d in 1:2) {
 z<-scenario_table[scenario_table$department==d,]
 plot(z$mean,1:nrow(z),xlim=range(c(z$lo90,z$hi90)),yaxt="n",ylab="",
  xlab="G1 minus G2 (percentage points)",main=paste("Department",d),pch=19)
 axis(2,at=1:nrow(z),labels=z$scenario,las=1,cex.axis=.8)
 segments(z$lo90,1:nrow(z),z$hi90,1:nrow(z),col="steelblue",lwd=2)
 abline(v=0,lty=2)
}
invisible(dev.off())
cat("END_19\n")

# Section 21
cat("BEGIN_21\n")
# Separate probability-scale teaching model for a transparent tipping point.
# Observed difference = standardized difference + h*(q1-q2).
observed_difference <- .06
latent_risk_difference <- .20
imbalance <- c(0,.1,.2,.3,.4)
print(data.frame(imbalance=imbalance,
  corrected_difference=observed_difference-latent_risk_difference*imbalance),
  row.names=FALSE)
cat("Sign tipping imbalance:",observed_difference/latent_risk_difference,"\n")
cat("For a 2-point decision threshold:",
    (observed_difference-.02)/latent_risk_difference,"\n")
cat("END_21\n")

# Section 25
cat("BEGIN_25\n")
# Real counts available in base R: admission x gender x department.
ucb <- as.data.frame(UCBAdmissions)
print(aggregate(Freq~Admit,ucb,sum),row.names=FALSE)
expanded <- ucb[rep(seq_len(nrow(ucb)),ucb$Freq),c("Admit","Gender","Dept")]
cat("Expanded rows:",nrow(expanded),"\n")
stopifnot(all(xtabs(~Admit+Gender+Dept,expanded)==UCBAdmissions))
cat("Columns after expansion:",paste(names(expanded),collapse=", "),"\n")
cat("END_25\n")

# Section 26
cat("BEGIN_26\n")
softmax <- function(x) {w<-exp(x-max(x));w/sum(w)}
base <- rep(0,6)
shift <- base;shift[1]<-1
print(data.frame(department=LETTERS[1:6],
  baseline=round(softmax(base),3),shifted=round(softmax(shift),3)),row.names=FALSE)
stopifnot(max(abs(softmax(shift)-softmax(shift+7)))<1e-12)
cat("END_26\n")

# Section 29
cat("BEGIN_29\n")
# Known loadings = 1, known independent noise SD = 1, prior u ~ N(0,1).
# Conjugate toy measurement model, not fitted UCB applicant measurements.
measurements <- c(1.2,.8,1.0)
proxy_update <- function(t) {
  variance<-1/(1+length(t)); c(mean=variance*sum(t),sd=sqrt(variance))
}
print(round(rbind(One_proxy=proxy_update(measurements[1]),
                  Three_proxies=proxy_update(measurements)),3))
cat("If all three are copies of one measurement, use one likelihood term.\n")
cat("END_29\n")

# Section 32
cat("BEGIN_32\n")
# X independent of binary U, with the same 50:50 mix in both X groups.
alpha <- -2; beta <- log(2); gamma <- 2
conditional <- outer(0:1,0:1,function(x,u) plogis(alpha+beta*x+gamma*u))
marginal <- rowMeans(conditional)
print(round(data.frame(X=0:1,U0=conditional[,1],U1=conditional[,2],
                       marginal=marginal),4),row.names=FALSE)
cat("Conditional OR:",exp(beta),"; marginal OR:",
    round(exp(diff(qlogis(marginal))),4),"\n")
cat("END_32\n")

# Section 34
cat("BEGIN_34\n")
# New applicants: draw fresh u, then D, then A using that simulated D.
# Parameters come from the fixed-scenario posterior, not from their priors.
set.seed(1034)
R <- 500
rep_summary <- matrix(NA_real_,R,4)
colnames(rep_summary)<-c("D2_G1","D2_G2","A_G1","A_G2")
for(s in 1:R) for(group in 1:2) {
 f<-fits$Official_fixed[[group]]
 n<-sum(G==group);u<-rnorm(n)
 d<-1+rbinom(n,1,plogis(f$delta[s]+f$g*u))
 a<-rbinom(n,1,plogis(f$a[s,d]+f$b*u))
 rep_summary[s,group]<-mean(d==2)
 rep_summary[s,group+2]<-mean(a)
}
observed<-c(mean(D[G==1]==2),mean(D[G==2]==2),mean(A[G==1]),mean(A[G==2]))
print(round(data.frame(observed=observed,t(apply(rep_summary,2,summary90))),3))
cat("END_34\n")

# Section 35
cat("BEGIN_35\n")
# One full prior draw for two groups; repeat before assessing plausibility.
set.seed(1035)
n <- 500
u_prior <- rnorm(n);group <- rep(1:2,each=n/2)
a_prior <- matrix(rnorm(4),2,2);delta_prior <- rnorm(2)
b_prior <- runif(2);g_prior <- runif(2)
d_prior <- 1+rbinom(n,1,plogis(delta_prior[group]+g_prior[group]*u_prior))
p_prior <- plogis(a_prior[cbind(group,d_prior)]+b_prior[group]*u_prior)
y_prior <- rbinom(n,1,p_prior)
print(data.frame(G=1:2,dept2=tapply(d_prior==2,group,mean),
                 admission=tapply(y_prior,group,mean)),row.names=FALSE)
cat("END_35\n")

# Section 36
cat("BEGIN_36\n")
# Artificial MCMC arrays; no claims about an HMC fit follow from these.
set.seed(1036)
x<-matrix(rnorm(4000),1000,4);bad<-x;bad[,4]<-bad[,4]+2
report<-function(z)c(Rhat=posterior::rhat(z),ESS_bulk=posterior::ess_bulk(z),
                    ESS_tail=posterior::ess_tail(z),MCSE=posterior::mcse_mean(z))
print(round(rbind(Agreeing=report(x),Shifted_chain=report(bad)),3))
cat("Posterior SD .03, mean-specific ESS 100: MCSE about",.03/sqrt(100),"\n")
cat("END_36\n")

# Section 43
cat("BEGIN_43\n")
# Invented business inputs, not causal estimates from a fitted campaign.
margin <- 40;cost <- .50
uplift <- c(.03,.015,.005)
print(data.frame(scenario=c("Optimistic","Intermediate","Pessimistic"),
  uplift=uplift,incremental_profit_per_customer=margin*uplift-cost),row.names=FALSE)
cat("Break-even conversion uplift:",cost/margin,"\n")
cat("END_43\n")

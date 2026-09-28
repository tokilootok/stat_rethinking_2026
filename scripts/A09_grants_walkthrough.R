# A09 grant laboratory, readable executable route. Run from repository root:
# Rscript scripts/A09_grants_walkthrough.R
# Requires rethinking only for NWOGrants. All fitting uses base-R grid integration.

# Section 53
cat("BEGIN_53\n")
# Start here. Run all blocks in order or run the complete script.
# rethinking is used only to load the dataset; fitting below uses base R.
data("NWOGrants",package="rethinking")
d <- data.frame(discipline_id=as.integer(NWOGrants$discipline),
  discipline=as.character(NWOGrants$discipline),
  gender_id=ifelse(NWOGrants$gender=="f",1L,2L),
  gender=ifelse(NWOGrants$gender=="f","female","male"),
  applications=NWOGrants$applications,awards=NWOGrants$awards)
d <- d[order(d$discipline_id,d$gender_id),]; rownames(d)<-NULL
keys <- paste(d$gender_id,d$discipline_id,sep="_")
ids <- sort(unique(d$discipline_id))
female <- match(paste(1,ids,sep="_"),keys)
male <- match(paste(2,ids,sep="_"),keys)
stopifnot(!anyNA(c(female,male)),!anyDuplicated(keys))
cat("Cells:",nrow(d),"; applications:",sum(d$applications),
    "; awards:",sum(d$awards),"\n")
cat("END_53\n")

# Section 55
cat("BEGIN_55\n")
totals <- aggregate(cbind(awards,applications)~gender,d,sum)
totals$rate <- totals$awards/totals$applications
print(totals,digits=4,row.names=FALSE)
cat("Raw female-minus-male gap (percentage points):",
    round(100*(totals$rate[totals$gender=="female"]-
               totals$rate[totals$gender=="male"]),3),"\n")
# Same number of awards, very different opportunities: a separate toy example.
print(data.frame(awards=c(10,10),applications=c(20,100),rate=c(.5,.1)))
cat("END_55\n")

# Section 57
cat("BEGIN_57\n")
# First actual cell: same sequence likelihood up to the binomial coefficient.
k<-1; n<-d$applications[k]; a<-d$awards[k]
y<-c(rep(1,a),rep(0,n-a)); probs<-c(.10,.25,.50)
comparison<-data.frame(p=probs,
  log_individual=sapply(probs,function(p) sum(dbinom(y,1,p,log=TRUE))),
  log_grouped=dbinom(a,n,probs,log=TRUE))
comparison$gap<-comparison$log_grouped-comparison$log_individual
print(round(comparison,4),row.names=FALSE)
stopifnot(max(abs(comparison$gap-lchoose(n,a)))<1e-10)
cat("Constant gap equals log choose(n,a):",round(lchoose(n,a),4),"\n")
cat("END_57\n")

# Section 59
cat("BEGIN_59\n")
# Teaching population with no confounding: X and U are independent.
p <- outer(0:1,0:1,function(x,u) plogis(-2+log(2)*x+2*u))
marginal <- rowMeans(p)
cat("Conditional odds ratio:",2,"\n")
cat("Marginal odds ratio:",round(exp(diff(qlogis(marginal))),3),"\n")
cat("Marginal probability difference:",round(diff(marginal),3),"\n")
cat("END_59\n")

# Section 60
cat("BEGIN_60\n")
rate <- d$awards/d$applications
nf <- d$applications[female]; nm <- d$applications[male]
w_f <- nf/sum(nf); w_m <- nm/sum(nm); w_pool <- (nf+nm)/sum(nf+nm)
names(w_f)<-names(w_m)<-names(w_pool)<-as.character(ids)
raw_total <- sum(w_f*rate[female])-sum(w_m*rate[male])
raw_standardized <- sum(w_pool*(rate[female]-rate[male]))
print(round(100*c(Observed_composition=raw_total,
                  Common_pooled_composition=raw_standardized),3))
cat("END_60\n")

# Section 61
cat("BEGIN_61\n")
weights <- data.frame(id=ids,discipline=d$discipline[female],
  female=w_f,male=w_m,pooled=w_pool,equal=rep(1/length(ids),length(ids)))
print(weights,digits=3,row.names=FALSE)
stopifnot(all(abs(colSums(weights[c("female","male","pooled","equal")])-1)<1e-12))
cat("END_61\n")

# Section 62
cat("BEGIN_62\n")
stopifnot(nrow(d)==18,length(ids)==9,sum(d$applications)==2823,
  sum(d$awards)==467,all(d$applications>0),all(d$awards>=0),
  all(d$awards<=d$applications),all(table(d$discipline_id)==2),
  !anyDuplicated(keys))
print(d,row.names=FALSE)
cat("PASS: 18 unique cells, nine complete discipline pairs, valid counts and totals.\n")
cat("END_62\n")

# Section 64
cat("BEGIN_64\n")
# Saturated binomial model: independent Normal(0,1) cell logits.
# Fit each one-dimensional posterior on an evenly spaced logit grid.
fit_cells <- function(y,n,prior_mean=0,prior_sd=1,step=.01,S=60000) {
  alpha<-seq(-10,10,by=step);p<-plogis(alpha)
  mass<-vapply(seq_along(y),function(k) {
    z<-dbinom(y[k],n[k],p,log=TRUE)+dnorm(alpha,prior_mean,prior_sd,log=TRUE)
    w<-exp(z-max(z));w/sum(w)
  },numeric(length(alpha)))
  stopifnot(max(colSums(mass[abs(alpha)>9,,drop=FALSE]))<1e-8)
  draws<-vapply(seq_along(y),function(k)
    sample(p,S,replace=TRUE,prob=mass[,k]),numeric(S))
  list(draws=draws,mean=colSums(mass*p))
}
summary90 <- function(x)c(mean=mean(x),lo90=unname(quantile(x,.05)),
  hi90=unname(quantile(x,.95)),Pr_positive=mean(x>0))
contrasts <- function(p) {
  pf<-p[,female,drop=FALSE];pm<-p[,male,drop=FALSE]
  cbind(Total=drop(pf%*%w_f-pm%*%w_m),
        Standardized=drop((pf-pm)%*%w_pool))
}
set.seed(9064)
fit_grid <- fit_cells(d$awards,d$applications)
p_draw <- fit_grid$draws
colnames(p_draw)<-keys
# Refine the grid independently of posterior sampling variation.
fine <- fit_cells(d$awards,d$applications,step=.005,S=10)
stopifnot(max(abs(fit_grid$mean-fine$mean))<1e-6)
cat("Fitted probabilities:",ncol(p_draw),"; posterior draws:",nrow(p_draw),"\n")
cat("Maximum cell-mean change after grid refinement:",
    format(max(abs(fit_grid$mean-fine$mean)),digits=3),"\n")
cat("END_64\n")

# Section 66
cat("BEGIN_66\n")
posterior_contrasts<-contrasts(p_draw)
posterior_table<-t(apply(posterior_contrasts,2,summary90))
print(round(posterior_table,4))
probabilities<-rbind(
  Observed_composition=c(Female=mean(p_draw[,female]%*%w_f),
                         Male=mean(p_draw[,male]%*%w_m)),
  Common_composition=c(Female=mean(p_draw[,female]%*%w_pool),
                       Male=mean(p_draw[,male]%*%w_pool)))
print(round(probabilities,4))
dir.create("courses/figures/A09_grants",recursive=TRUE,showWarnings=FALSE)
write.csv(posterior_table,"courses/figures/A09_grants/posterior_summary.csv")
png("courses/figures/A09_grants/contrasts.png",width=1100,height=550,res=130)
z<-100*posterior_table
par(mar=c(5,8,3,1))
plot(z[,"mean"],1:2,xlim=range(z[,c("lo90","hi90")]),ylim=c(.6,2.4),
  yaxt="n",ylab="",xlab="Female minus male (percentage points)",pch=19,
  main="NWOGrants: changing the target composition")
axis(2,at=1:2,labels=rownames(z),las=1)
segments(z[,"lo90"],1:2,z[,"hi90"],1:2,col="steelblue",lwd=3)
points(z[,"mean"],1:2,pch=19);abline(v=0,lty=2)
invisible(dev.off())
cat("END_66\n")

# Section 68
cat("BEGIN_68\n")
# One transparent decomposition using the common pooled composition.
pf<-fit_grid$mean[female];pm<-fit_grid$mean[male]
common<-sum(w_pool*(pf-pm))
female_shift<-sum((w_f-w_pool)*pf)
male_shift<- -sum((w_m-w_pool)*pm)
reconstructed<-common+female_shift+male_shift
print(round(100*c(Common_composition=common,
  Female_composition_shift=female_shift,Male_composition_shift=male_shift,
  Total_reconstructed=reconstructed),3))
stopifnot(abs(reconstructed-(sum(w_f*pf)-sum(w_m*pm)))<1e-12)
cat("END_68\n")

# Section 69
cat("BEGIN_69\n")
correct<-sum(w_f*pf)-sum(w_m*pm)
# Deliberately attach both gender-specific weight vectors in reverse ID order.
wrong<-sum(rev(w_f)*pf)-sum(rev(w_m)*pm)
print(round(100*c(Correct=correct,Permuted_weights=wrong),3))
# A safe function demands matching keys instead of trusting vector position.
safe_average<-function(p,w) {
  stopifnot(!is.null(names(p)),!is.null(names(w)),
            setequal(names(p),names(w)))
  sum(p*w[names(p)])
}
names(pf)<-names(pm)<-as.character(ids)
recovered<-safe_average(pf,rev(w_f))-safe_average(pm,rev(w_m))
stopifnot(abs(correct-recovered)<1e-12)
cat("PASS: identifier alignment recovers the result even when weight order changes.\n")
cat("END_69\n")

# Section 70
cat("BEGIN_70\n")
set.seed(9070)
prior_counts<-function(mu,sd) {
  pp<-matrix(plogis(rnorm(10000*nrow(d),mu,sd)),ncol=nrow(d))
  yy<-matrix(rbinom(length(pp),size=rep(d$applications,each=nrow(pp)),
                    prob=as.vector(pp)),nrow=nrow(pp))
  rowSums(yy)
}
prior_default<-prior_counts(0,1)
prior_centered<-prior_counts(qlogis(.17),.75)
print(round(rbind(Normal_0_1=summary90(prior_default)[1:3],
  Rate_centered=summary90(prior_centered)[1:3]),1))
cat("Observed total awards:",sum(d$awards),"\n")
cat("END_70\n")

# Section 71
cat("BEGIN_71\n")
# One cell shows expected-count uncertainty versus new count variation.
set.seed(9071)
k<-which.min(d$applications)
expected<-d$applications[k]*p_draw[,k]
replicated<-rbinom(nrow(p_draw),d$applications[k],p_draw[,k])
print(d[k,c("discipline","gender","applications","awards")],row.names=FALSE)
print(round(rbind(Expected_count=summary90(expected)[1:3],
                  Replicated_count=summary90(replicated)[1:3]),3))
cat("END_71\n")

# Section 72
cat("BEGIN_72\n")
set.seed(9072)
y_rep<-vapply(seq_len(nrow(d)),function(k)
  rbinom(nrow(p_draw),d$applications[k],p_draw[,k]),numeric(nrow(p_draw)))
colnames(y_rep)<-keys
rep_totals<-cbind(All=rowSums(y_rep),Female=rowSums(y_rep[,female]),
                  Male=rowSums(y_rep[,male]))
observed<-c(All=sum(d$awards),Female=sum(d$awards[female]),Male=sum(d$awards[male]))
print(round(data.frame(observed,t(apply(rep_totals,2,function(x)
  summary90(x)[1:3]))),1))
# Pearson discrepancies use the same posterior p for observed and replicated data.
expected<-sweep(p_draw,2,d$applications,"*")
variance<-sweep(p_draw*(1-p_draw),2,d$applications,"*")
observed_matrix<-matrix(d$awards,nrow(p_draw),nrow(d),byrow=TRUE)
T_obs<-rowSums((observed_matrix-expected)^2/variance)
T_rep<-rowSums((y_rep-expected)^2/variance)
cat("Posterior predictive tail fraction for Pearson discrepancy:",
    round(mean(T_rep>=T_obs),3),"\n")
png("courses/figures/A09_grants/ppc.png",width=1000,height=600,res=130)
hist(rep_totals[,"All"],breaks=35,col="lightblue",border="white",
  main="NWOGrants: replicated total awards",xlab="Total awards in a replicated cohort")
abline(v=sum(d$awards),col="firebrick",lwd=3)
legend("topright","Observed total",col="firebrick",lwd=3,bty="n")
invisible(dev.off())
cat("END_72\n")

# Section 73
cat("BEGIN_73\n")
# One recovery exercise: no within-discipline difference, unequal composition.
# Chosen truths, not prior draws: this is NOT SBC.
set.seed(9073)
base<-seq(.08,.30,length.out=length(ids))
p_true<-base[match(d$discipline_id,ids)]
y_sim<-rbinom(nrow(d),d$applications,p_true)
recovery<-fit_cells(y_sim,d$applications,S=30000)
true_contrast<-contrasts(matrix(p_true,nrow=1))
recovered_contrast<-contrasts(recovery$draws)
print(round(cbind(truth=as.numeric(true_contrast),
  t(apply(recovered_contrast,2,summary90))[,1:3]),4))
stopifnot(abs(true_contrast[1,"Standardized"])<1e-12)
cat("One recovery run checks the pipeline; it does not establish coverage or SBC.\n")
cat("END_73\n")

# Section 75
cat("BEGIN_75\n")
set.seed(9075)
alternative<-fit_cells(d$awards,d$applications,prior_mean=qlogis(.17),prior_sd=.75)
scenarios<-list(Default_total=posterior_contrasts[,"Total"],
  Default_pooled=posterior_contrasts[,"Standardized"],
  Default_equal=rowMeans(p_draw[,female]-p_draw[,male]),
  Centered_prior_pooled=contrasts(alternative$draws)[,"Standardized"])
print(round(t(vapply(scenarios,summary90,numeric(4))),4))
# Same independent cell fit, excluding one discipline and renormalizing weights.
# This changes the target population; it is not new-discipline prediction.
leave_out<-sapply(seq_along(ids),function(j) {
  w<-w_pool[-j]/sum(w_pool[-j])
  mean((p_draw[,female[-j],drop=FALSE]-p_draw[,male[-j],drop=FALSE])%*%w)
})
cat("Leave-one-discipline-out standardized posterior means, range:",
    paste(round(range(leave_out),4),collapse=" to "),"\n")
cat("END_75\n")

# Section 76
cat("BEGIN_76\n")
# Illustration only: a shared random cell probability creates extra variation.
n<-100;mu<-.17;kappa<-20
alpha<-mu*kappa;beta<-(1-mu)*kappa
binomial_variance<-n*mu*(1-mu)
beta_binomial_variance<-binomial_variance*(n+kappa)/(1+kappa)
print(data.frame(model=c("Binomial","Beta-binomial"),mean=n*mu,
                 variance=c(binomial_variance,beta_binomial_variance)),digits=4,
      row.names=FALSE)
cat("END_76\n")

# Section 78
cat("BEGIN_78\n")
checks<-c(complete_cells=nrow(d)==18 && !anyDuplicated(keys),
  valid_counts=all(d$awards>=0 & d$awards<=d$applications),
  correct_totals=sum(d$applications)==2823 && sum(d$awards)==467,
  weights_normalized=all(abs(c(sum(w_f),sum(w_m),sum(w_pool))-1)<1e-12),
  aligned_draws=identical(colnames(p_draw),keys),
  refinement=max(abs(fit_grid$mean-fine$mean))<1e-6,
  decomposition=abs(reconstructed-correct)<1e-12,
  permutation_defense=abs(correct-recovered)<1e-12)
print(checks)
stopifnot(all(checks))
cat("PASS for these deterministic checks only; HMC and cross-language parity not run.\n")
cat("END_78\n")

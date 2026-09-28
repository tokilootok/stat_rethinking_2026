# A09 Problem 1A walkthrough. Run from the repository root.
# Rscript scripts/A09_problem1a_walkthrough.R
# Shared numerical implementation: scripts/A09_walkthrough_helpers.R

# Lesson section 4: Generate probabilities and binary observations
cat("BEGIN_4\n")
source("scripts/A09_walkthrough_helpers.R")
set.seed(9104)
x<-runif(120,-2,2)
alpha_true<- -1;beta_true<-1
p_true<-plogis(alpha_true+beta_true*x)
y<-rbinom(length(x),1,p_true)
print(round(data.frame(x=x[1:6],eta=alpha_true+beta_true*x[1:6],
  probability=p_true[1:6],outcome=y[1:6]),3),row.names=FALSE)
cat("Events:",sum(y),"out of",length(y),"\n")
cat("END_4\n")

# Lesson section 7: Check individual versus grouped binary likelihoods
cat("BEGIN_7\n")
# Aggregation is exact here only within a common probability, not across x values.
z<-c(1,0,0,1,0);p<-.3
print(c(log_individual=sum(dbinom(z,1,p,log=TRUE)),
        log_grouped=dbinom(sum(z),length(z),p,log=TRUE),
        difference=lchoose(length(z),sum(z))))
stopifnot(abs(dbinom(sum(z),length(z),p,log=TRUE)-
  sum(dbinom(z,1,p,log=TRUE))-lchoose(length(z),sum(z)))<1e-12)
cat("END_7\n")

# Lesson section 8: Simulate the probability implications of intercept priors
cat("BEGIN_8\n")
set.seed(9108)
prior_moments<-t(sapply(c(.5,1,2),function(sd) {
  probability<-plogis(rnorm(50000,0,sd))
  c(mean=mean(probability),variance=var(probability),
    lo90=unname(quantile(probability,.05)),hi90=unname(quantile(probability,.95)))
}))
rownames(prior_moments)<-c("alpha_SD_0.5","alpha_SD_1","alpha_SD_2")
print(round(prior_moments,4))
cat("END_8\n")

# Lesson section 9: Combine intercept and slope priors before reading an effect
cat("BEGIN_9\n")
set.seed(9109)
a<-rnorm(50000,0,1);b<-rnorm(50000,0,1)
prior_difference<-plogis(a+b)-plogis(a)
print(round(rbind(Probability_change_0_to_1=a09_summary(prior_difference),
                  Odds_multiplier=a09_summary(exp(b))),3))
cat("END_9\n")

# Lesson section 10: Update the same sample sequentially
cat("BEGIN_10\n")
set.seed(9110)
seen<-c(0,10,40,120)
fits<-lapply(seen,function(n) a09_regression(x[seq_len(n)],y[seq_len(n)]))
posterior<-lapply(fits,a09_draw)
learning<-do.call(rbind,lapply(seq_along(seen),function(j)
  data.frame(observed=seen[j],parameter=c("alpha","beta"),
    t(vapply(posterior[[j]],a09_summary,numeric(3))))))
print(round_table<-transform(learning,mean=round(mean,3),lo90=round(lo90,3),hi90=round(hi90,3)),
      row.names=FALSE)
post<-posterior[[length(seen)]]
fine<-a09_regression(x,y,agrid=seq(-7,7,length.out=401),bgrid=seq(-6,6,length.out=401))
mean_coef<-function(f) colSums(f$parameters*f$mass)
stopifnot(max(abs(mean_coef(fits[[4]])-mean_coef(fine)))<.001)
cat("PASS: posterior means stable within 0.001 after doubling grid resolution.\n")
cat("END_10\n")

# Lesson section 12: Read posterior probability contrasts and predictions
cat("BEGIN_12\n")
p0<-plogis(post$alpha);p1<-plogis(post$alpha+post$beta)
contrasts<-cbind(Probability_difference=p1-p0,Probability_ratio=p1/p0,
                 Odds_ratio=exp(post$beta))
print(round(t(apply(contrasts,2,a09_summary)),3))
set.seed(9112)
future<-rbinom(nrow(post),100,p1)
print(round(rbind(Expected_events=a09_summary(100*p1),
                  Future_events=a09_summary(future)),2))
out<-a09_output("1a")
write.csv(learning,file.path(out,"learning.csv"),row.names=FALSE)
png(file.path(out,"learning.png"),width=1200,height=600,res=130)
par(mfrow=c(1,2),mar=c(5,4,3,1))
xx<-seq(-2,2,length.out=81)
curve<-sapply(xx,function(v) plogis(post$alpha+post$beta*v))
band<-apply(curve,2,quantile,c(.05,.95))
plot(xx,colMeans(curve),type="l",ylim=c(0,1),xlab="Predictor x",ylab="Event probability",
 main="Posterior curve and simulation truth",lwd=2)
polygon(c(xx,rev(xx)),c(band[1,],rev(band[2,])),col=adjustcolor("steelblue",.2),border=NA)
lines(xx,colMeans(curve),lwd=2);lines(xx,plogis(alpha_true+beta_true*xx),col="firebrick",lty=2)
legend("topleft",c("Posterior mean","Generating probability"),lty=c(1,2),
 col=c("black","firebrick"),bty="n")
hist(100*(p1-p0),breaks=30,col="lightblue",border="white",
 main="Change from x = 0 to x = 1",xlab="Probability difference (percentage points)")
abline(v=0,lty=2);invisible(dev.off())
cat("END_12\n")

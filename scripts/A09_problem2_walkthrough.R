# A09 Problem 2 walkthrough. Run from the repository root.
# Rscript scripts/A09_problem2_walkthrough.R
# Shared numerical implementation: scripts/A09_walkthrough_helpers.R

# Lesson section 25: Translate a log link and simulate prior expected counts
cat("BEGIN_25\n")
source("scripts/A09_walkthrough_helpers.R")
set.seed(9325)
a<-rnorm(30000,3,.5);b<-rnorm(30000,0,.2)
prior_table<-t(sapply(c(-1,0,1),function(x) a09_summary(exp(a+b*x))))
rownames(prior_table)<-c("x_minus_1","x_zero","x_plus_1")
print(round(prior_table,2))
cat("A slope log(1.2) multiplies the expected count by",exp(log(1.2)),"per predictor unit.\n")
cat("END_25\n")

# Lesson section 29: Inspect the actual Kline observation unit
cat("BEGIN_29\n")
data("Kline",package="rethinking")
d<-Kline
d$x<-as.numeric(scale(log(d$population)))
print(d[c("culture","population","contact","total_tools","x")],digits=3,row.names=FALSE)
cat("Observations:",nrow(d),"societies; tools:",sum(d$total_tools),"\n")
cat("END_29\n")

# Lesson section 30: Fit four count models and compare numerical leave-one-out predictions
cat("BEGIN_30\n")
# All models use alpha ~ Normal(3,.5), slopes ~ Normal(0,.2).
fit_poisson<-function(rows,slope=TRUE) a09_regression(d$x[rows],d$total_tools[rows],
 family="poisson",prior_mean=c(3,0),prior_sd=c(.5,.2),
 agrid=seq(0,6,length.out=241),bgrid=if(slope) seq(-1.5,1.5,length.out=241) else 0)
all<-seq_len(nrow(d));groups<-split(all,d$contact)
models<-list(Intercept=list(rows=list(all),fits=list(fit_poisson(all,FALSE))),
 Population=list(rows=list(all),fits=list(fit_poisson(all))),
 Contact=list(rows=groups,fits=lapply(groups,fit_poisson,slope=FALSE)),
 Interaction=list(rows=groups,fits=lapply(groups,fit_poisson)))
loo<-sapply(models,function(m) {
 z<-numeric(nrow(d));for(j in seq_along(m$fits)) z[m$rows[[j]]]<-a09_loo(m$fits[[j]])
 z
})
difference<-sweep(loo,1,loo[,"Population"],"-")
print(data.frame(model=colnames(loo),elpd_loo=round(colSums(loo),2),
  difference_from_population=round(colSums(difference),2),
  paired_SE=round(sqrt(nrow(d)*apply(difference,2,var)),2)),row.names=FALSE)
print(data.frame(culture=d$culture,population_loo=round(loo[,"Population"],2),
  interaction_loo=round(loo[,"Interaction"],2)),row.names=FALSE)
cat("Higher elpd is better for leaving out one society; these are quadrature LOO scores, not PSIS.\n")
set.seed(9330)
post<-lapply(models$Interaction$fits,a09_draw)
coef_table<-do.call(rbind,lapply(names(post),function(g)
 data.frame(contact=g,parameter=c("alpha","beta"),
  round(t(vapply(post[[g]],a09_summary,numeric(3))),3))))
print(coef_table,row.names=FALSE)
# Resolution check on both interaction groups.
for(g in names(groups)) {
 fine<-a09_regression(d$x[groups[[g]]],d$total_tools[groups[[g]]],family="poisson",
  prior_mean=c(3,0),prior_sd=c(.5,.2),agrid=seq(0,6,length.out=481),
  bgrid=seq(-1.5,1.5,length.out=481))
 coarse<-models$Interaction$fits[[g]]
 stopifnot(max(abs(colSums(fine$parameters*fine$mass)-
                   colSums(coarse$parameters*coarse$mass)))<.001)
}
out<-a09_output("2")
write.csv(data.frame(culture=d$culture,loo),file.path(out,"loo_by_society.csv"),row.names=FALSE)
cat("END_30\n")

# Lesson section 34: Distinguish expected tool counts from replicated tool counts
cat("BEGIN_34\n")
# Section 34 is inside the current combined 30–35 lesson block.
set.seed(9334)
population<-10000
xx<-(log(population)-mean(log(d$population)))/sd(log(d$population))
forecast<-do.call(rbind,lapply(names(post),function(g) {
 lambda<-exp(post[[g]]$alpha+post[[g]]$beta*xx)
 data.frame(contact=g,quantity=c("Expected tools","Replicated tools"),
  round(rbind(a09_summary(lambda),a09_summary(rpois(length(lambda),lambda))),2))
}))
print(forecast,row.names=FALSE)
# Posterior predictive replication for all observed societies.
y_rep<-sapply(seq_len(nrow(d)),function(i) {
 pp<-post[[as.character(d$contact[i])]]
 rpois(nrow(pp),exp(pp$alpha+pp$beta*d$x[i]))
})
print(data.frame(culture=d$culture,observed=d$total_tools,
 round(t(apply(y_rep,2,a09_summary)),1)),row.names=FALSE)
png(file.path(out,"predictions.png"),width=1200,height=650,res=130)
par(mfrow=c(1,2),mar=c(5,4,3,1))
for(g in names(groups)) {
 r<-groups[[g]];xx<-seq(min(d$x[r]),max(d$x[r]),length.out=80)
 curves<-sapply(xx,function(z) exp(post[[g]]$alpha+post[[g]]$beta*z))
 band<-apply(curves,2,quantile,c(.05,.95))
 raw_pop<-exp(xx*sd(log(d$population))+mean(log(d$population)))
 plot(raw_pop,colMeans(curves),type="n",log="x",ylim=range(c(band,d$total_tools[r])),
  xlab="Population (log axis)",ylab="Expected tool count",main=paste(g,"contact"))
 polygon(c(raw_pop,rev(raw_pop)),c(band[1,],rev(band[2,])),col=adjustcolor("steelblue",.2),border=NA)
 lines(raw_pop,colMeans(curves),lwd=2);points(d$population[r],d$total_tools[r],pch=19)
}
invisible(dev.off())
cat("END_34\n")

# Lesson section 36: Read the innovation–loss model without claiming a new fit
cat("BEGIN_36\n")
# Hypothetical process parameters: expected equilibrium tools = alpha * P^gamma / loss.
P<-c(1000,10000,100000);innovation<-.5;gamma<-.5;loss<-1
print(data.frame(population=P,expected_tools=innovation*P^gamma/loss,
 doubled_innovation=2*innovation*P^gamma/loss,
 doubled_loss=innovation*P^gamma/(2*loss)),digits=3,row.names=FALSE)
cat("Only the innovation/loss ratio is identified from this mean curve without more assumptions.\n")
cat("END_36\n")

# Lesson section 37: Quantify the rare-event approximation
cat("BEGIN_37\n")
lambda<-4
print(data.frame(trials=c(20,100,1000),probability=lambda/c(20,100,1000),
 binomial_zero=dbinom(0,c(20,100,1000),lambda/c(20,100,1000)),
 poisson_zero=dpois(0,lambda)),digits=5,row.names=FALSE)
cat("END_37\n")

# Lesson section 38: Compare counts with unequal exposures
cat("BEGIN_38\n")
events<-c(6,3);days<-c(30,10)
print(data.frame(group=c("A","B"),events,days,rate=events/days,
 expected_at_30_days=30*events/days),row.names=FALSE)
cat("Offset: log(days) has coefficient one; it is not an estimated slope.\n")
cat("END_38\n")

# Lesson section 39: Make extra variation explicit
cat("BEGIN_39\n")
mu<-4;phi<-2
print(data.frame(model=c("Poisson","Gamma-Poisson"),mean=mu,
 variance=c(mu,mu+mu^2/phi),zero_probability=c(dpois(0,mu),dnbinom(0,size=phi,mu=mu))),
 digits=4,row.names=FALSE)
cat("A variance larger than the pooled mean is not alone a conditional model failure.\n")
cat("END_39\n")

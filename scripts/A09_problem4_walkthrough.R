# A09 Problem 4 walkthrough. Run from the repository root.
# Rscript scripts/A09_problem4_walkthrough.R
# Shared numerical implementation: scripts/A09_walkthrough_helpers.R

# Lesson section 42a: Generate one coherent campaign and its repeated purchases
cat("BEGIN_42a\n")
source("scripts/A09_walkthrough_helpers.R")
set.seed(9442)
N<-1200;A<-rbinom(N,1,.5);segment<-rbinom(N,1,.5)
rate_true<-.01*exp(.35*A+.4*segment)
# Independent increments from the same Poisson purchase process.
C15<-rpois(N,15*rate_true);C30<-C15+rpois(N,15*rate_true)
C60<-C30+rpois(N,30*rate_true);Y30<-as.integer(C30>0)
# A second observation scheme: noninformative follow-up conditional on assignment.
days<-ifelse(A==1,sample(c(15,30),N,replace=TRUE),sample(c(30,60),N,replace=TRUE))
observed_count<-ifelse(days==15,C15,ifelse(days==30,C30,C60))
customer<-data.frame(A,segment,Y30,C30,days,observed_count)
cells<-aggregate(cbind(customers=rep(1,N),buyers=Y30,purchases=C30,
                       observed_days=days,observed_purchases=observed_count)~A+segment,
                 customer,sum)
cells<-cells[order(cells$segment,cells$A),]
print(cells,row.names=FALSE)
i0<-which(cells$A==0);i1<-which(cells$A==1)
w<-as.numeric(table(factor(segment,levels=0:1)))/N
cat("Target segment weights:",w,"\n")
stopifnot(all(Y30<=C30),all(cells$buyers<=cells$customers))
cat("END_42a\n")

# Lesson section 42b: Estimate standardized conversion uplift
cat("BEGIN_42b\n")
set.seed(9443)
# Four independent cell logits; allows an assignment-by-segment interaction.
conversion<-a09_cells(cells$buyers,cells$customers,prior_mean=qlogis(.3),prior_sd=1)
conversion_uplift<-drop((conversion$draws[,i1]-conversion$draws[,i0])%*%w)
cat("Standardized conversion uplift (percentage points):\n")
print(round(100*a09_summary(conversion_uplift),2))
cat("Generating mean uplift for these customer segments (percentage points):",
 round(100*mean(exp(-30*.01*exp(.4*segment))-exp(-30*.01*exp(.35+.4*segment))),2),"\n")
cat("END_42b\n")

# Lesson section 43: Estimate purchases rather than buyers
cat("BEGIN_43\n")
set.seed(9444)
# Sum of independent Poisson counts: exposure here is number of 30-day customers.
count_fit<-a09_cells(cells$purchases,cells$customers,family="poisson",
                    prior_mean=log(.4),prior_sd=1)
count_uplift<-drop((count_fit$draws[,i1]-count_fit$draws[,i0])%*%w)
print(round(rbind(Extra_buyers_per_1000=a09_summary(1000*conversion_uplift),
                  Extra_purchases_per_1000=a09_summary(1000*count_uplift)),1))
cat("These are separate summaries of the same process, not independent datasets to multiply together.\n")
cat("END_43\n")

# Lesson section 44: Correct unequal follow-up with exposure
cat("BEGIN_44\n")
set.seed(9445)
# Now the modeled cell parameter is purchases per day; aggregate exposure is total days.
rate_fit<-a09_cells(cells$observed_purchases,cells$observed_days,family="poisson",
                   prior_mean=log(.01),prior_sd=1)
rate_uplift<-drop((rate_fit$draws[,i1]-rate_fit$draws[,i0])%*%w)
raw<-aggregate(cbind(observed_count,days)~A,customer,mean)
print(raw,digits=4,row.names=FALSE)
cat("Extra expected purchases per 1000 customers at common 30-day exposure:\n")
print(round(a09_summary(1000*30*rate_uplift),1))
cat("END_44\n")

# Lesson section 45: Verify what aggregation preserves
cat("BEGIN_45\n")
# Within an assignment/segment cell, the daily rate is common but exposure can differ.
r<-which(A==0 & segment==0)
candidate_rates<-c(.005,.01,.02)
individual<-sapply(candidate_rates,function(rate)
 sum(dpois(observed_count[r],days[r]*rate,log=TRUE)))
grouped<-dpois(sum(observed_count[r]),sum(days[r])*candidate_rates,log=TRUE)
print(round(data.frame(rate=candidate_rates,log_individual=individual,
 log_grouped=grouped,gap=grouped-individual),4),row.names=FALSE)
stopifnot(max(abs((grouped-individual)-(grouped-individual)[1]))<1e-8)
cat("Aggregation preserves this rate posterior up to a constant; it does not recover missing customer predictors.\n")
cat("END_45\n")

# Lesson section 46: Turn the defined count contrast into a transparent decision report
cat("BEGIN_46\n")
# Fixed illustrative margin and campaign cost, not estimated financial parameters.
margin_per_purchase<-20;cost_per_customer<-1.2
profit<-margin_per_purchase*count_uplift-cost_per_customer
report<-rbind(Conversion_difference=a09_summary(conversion_uplift),
 Purchases_per_customer_30d=a09_summary(count_uplift),
 Profit_per_customer=a09_summary(profit))
print(round(report,3))
cat("Pr(incremental profit > 0) under these fixed economic inputs:",round(mean(profit>0),3),"\n")
# Simple posterior predictive check at the original common 30-day exposure.
set.seed(9446)
rep_counts<-sapply(seq_len(nrow(cells)),function(k)
 rpois(nrow(count_fit$draws),cells$customers[k]*count_fit$draws[,k]))
print(round(c(observed_purchases=sum(C30),a09_summary(rowSums(rep_counts))),1))
out<-a09_output("4");write.csv(report,file.path(out,"campaign_report.csv"))
z<-rbind(Buyers=1000*a09_summary(conversion_uplift),Purchases=1000*a09_summary(count_uplift))
png(file.path(out,"campaign.png"),width=1000,height=550,res=130)
par(mar=c(5,7,3,1))
plot(z[,1],1:2,xlim=range(c(0,z[,2:3])),ylim=c(.6,2.4),yaxt="n",ylab="",pch=19,
 xlab="Additional expected outcomes per 1,000 customers / 30 days",main="One campaign, two outcomes")
axis(2,at=1:2,labels=rownames(z),las=1)
segments(z[,2],1:2,z[,3],1:2,col="steelblue",lwd=3);points(z[,1],1:2,pch=19)
abline(v=0,lty=2);invisible(dev.off())
cat("END_46\n")

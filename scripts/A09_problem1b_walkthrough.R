# A09 Problem 1B walkthrough. Run from the repository root.
# Rscript scripts/A09_problem1b_walkthrough.R
# Shared numerical implementation: scripts/A09_walkthrough_helpers.R

# Lesson section 14: Distinguish the small teaching example from real admissions
cat("BEGIN_14\n")
source("scripts/A09_walkthrough_helpers.R")
# Invented two-department example, retaining category 1 minus category 2.
print(data.frame(comparison=c("Overall","Within either department","Common 50:50 mix"),
  difference=c(.32-.58,.20-.10,.5*(.20+.80)-.5*(.10+.70))),row.names=FALSE)
# Real UCB counts in base R; now category 1 = male, category 2 = female.
tab<-as.data.frame(UCBAdmissions)
cells<-aggregate(Freq~Gender+Dept,tab,sum);names(cells)[3]<-"applications"
admitted<-tab[tab$Admit=="Admitted",]
cells$admitted<-admitted$Freq[match(paste(cells$Gender,cells$Dept),paste(admitted$Gender,admitted$Dept))]
cells$gid<-ifelse(cells$Gender=="Male",1,2)
cells<-cells[order(cells$Dept,cells$gid),]
keys<-paste(cells$gid,cells$Dept,sep="_")
ids<-levels(cells$Dept)
i1<-match(paste(1,ids,sep="_"),keys);i2<-match(paste(2,ids,sep="_"),keys)
totals<-aggregate(cbind(admitted,applications)~gid,cells,sum)
totals$rate<-totals$admitted/totals$applications
print(totals,digits=4,row.names=FALSE)
set.seed(9214)
aggregate_fit<-a09_cells(totals$admitted,totals$applications,prior_sd=1.5)
aggregate_difference<-aggregate_fit$draws[,1]-aggregate_fit$draws[,2]
cat("Aggregate posterior (male minus female), percentage points:\n")
print(round(100*a09_summary(aggregate_difference),2))
cat("END_14\n")

# Lesson section 15: Fit the twelve observed gender-department cells
cat("BEGIN_15\n")
set.seed(9215)
f<-a09_cells(cells$admitted,cells$applications,prior_sd=1.5)
p_draw<-f$draws
within<-p_draw[,i1]-p_draw[,i2];colnames(within)<-ids
print(round(100*t(apply(within,2,a09_summary)),2))
finer<-a09_cells(cells$admitted,cells$applications,prior_sd=1.5,step=.005,S=10)
stopifnot(max(abs(f$mean-finer$mean))<1e-6)
cat("END_15\n")

# Lesson section 18: Standardize both categories to the same department population
cat("BEGIN_18\n")
n1<-cells$applications[i1];n2<-cells$applications[i2]
w1<-n1/sum(n1);w2<-n2/sum(n2);w<-(n1+n2)/sum(n1+n2)
print(data.frame(department=ids,male=round(w1,3),female=round(w2,3),pooled=round(w,3)),
      row.names=FALSE)
comparisons<-cbind(Aggregate_model=aggregate_difference,
  Cell_model_own_mix=drop(p_draw[,i1]%*%w1-p_draw[,i2]%*%w2),
  Cell_model_pooled_mix=drop(within%*%w),Equal_departments=rowMeans(within))
print(round(100*t(apply(comparisons,2,a09_summary)),2))
# Raw own-mixture arithmetic must reproduce the observed aggregate rates.
stopifnot(abs(sum(w1*cells$admitted[i1]/n1)-totals$rate[1])<1e-12)
out<-a09_output("1b")
write.csv(t(apply(comparisons,2,a09_summary)),file.path(out,"ucb_contrasts.csv"))
z<-100*t(apply(comparisons,2,a09_summary))
png(file.path(out,"ucb_contrasts.png"),width=1150,height=600,res=130)
par(mar=c(5,11,3,1))
plot(z[,1],1:4,xlim=range(z[,2:3]),yaxt="n",ylab="",pch=19,
 main="Real admissions: different target comparisons",xlab="Male minus female (percentage points)")
axis(2,at=1:4,labels=rownames(z),las=1,cex.axis=.8)
segments(z[,2],1:4,z[,3],1:4,col="steelblue",lwd=3);points(z[,1],1:4,pch=19)
abline(v=0,lty=2);invisible(dev.off())
cat("END_18\n")

# Lesson section 22: Check predictions while preserving cell labels and denominators
cat("BEGIN_22\n")
set.seed(9222)
y_rep<-vapply(seq_len(nrow(cells)),function(k)
 rbinom(nrow(p_draw),cells$applications[k],p_draw[,k]),numeric(nrow(p_draw)))
summary<-cbind(All=rowSums(y_rep),Male=rowSums(y_rep[,i1]),Female=rowSums(y_rep[,i2]))
observed<-c(sum(cells$admitted),sum(cells$admitted[i1]),sum(cells$admitted[i2]))
print(round(data.frame(observed,t(apply(summary,2,a09_summary))),1))
print(data.frame(department=cells$Dept,gender=cells$Gender,observed=cells$admitted,
  round(t(apply(y_rep,2,a09_summary)),1)),row.names=FALSE)
cat("No causal effect or within-cell calibration is established by this predictive check.\n")
cat("END_22\n")

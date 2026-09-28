# A09 Problem 3 walkthrough. Run from the repository root.
# Rscript scripts/A09_problem3_walkthrough.R
# Shared numerical implementation: scripts/A09_walkthrough_helpers.R

# Lesson section 40: Separate non-collapsibility from confounding
cat("BEGIN_40\n")
source("scripts/A09_walkthrough_helpers.R")
# Exact population calculation: X and U independent, Pr(U=1)=0.5 in both X groups.
p<-outer(0:1,0:1,function(x,u) plogis(-2+log(2)*x+2*u))
marginal<-rowMeans(p)
print(round(data.frame(X=0:1,U0=p[,1],U1=p[,2],marginal=marginal),4),row.names=FALSE)
cat("Conditional OR:",2,"; marginal OR:",round(exp(diff(qlogis(marginal))),3),"\n")
cat("No unequal latent composition is present in this example.\n")
cat("END_40\n")

# Lesson section 41a: Compare observed association with a known intervention effect
cat("BEGIN_41a\n")
# Declared toy causal process: U -> treatment and U -> outcome, plus treatment -> outcome.
# Pr(U=1)=.5; Pr(T=1|U=0)=.2; Pr(T=1|U=1)=.8.
# Pr(Y=1|T,U)=.10+.10*T+.40*U, valid for all four combinations.
selected_u<-c(Control=.2,Treated=.8)
observed<-c(Control=.10+.40*selected_u[1],Treated=.20+.40*selected_u[2])
intervened<-c(Control=.10+.40*.5,Treated=.20+.40*.5)
print(data.frame(scenario=c("Control","Treated"),observed=as.numeric(observed),
 intervened=as.numeric(intervened)),row.names=FALSE)
association<-unname(diff(observed));effect<-unname(diff(intervened))
cat("Observed difference:",association,"; true intervention difference:",effect,"\n")
stopifnot(abs(association-.34)<1e-12,abs(effect-.10)<1e-12)
cat("END_41a\n")

# Lesson section 41b: Check why conditioning on a collider is different
cat("BEGIN_41b\n")
# Admissions mechanism used in A10: category 1 selects D2 only when U=1;
# category 2 selects D2 independently of U; original Pr(U=1)=.1 in both.
observed_d2<-c(G1=.30,G2=.9*.30+.1*.50)
controlled_d2<-c(G1=.9*.10+.1*.30,G2=.9*.30+.1*.50)
print(data.frame(comparison=c("Observed within D2","Intervene on D2, common U mix"),
 difference=c(observed_d2[1]-observed_d2[2],controlled_d2[1]-controlled_d2[2])),row.names=FALSE)
out<-a09_output("3")
comparisons<-data.frame(case=c("Confounded treatment","Selected admissions"),
 observed=c(association,unname(observed_d2[1]-observed_d2[2])),
 intervention=c(effect,unname(controlled_d2[1]-controlled_d2[2])))
write.csv(comparisons,file.path(out,"known_population_contrasts.csv"),row.names=FALSE)
png(file.path(out,"causal_contrasts.png"),width=1000,height=550,res=130)
par(mar=c(5,10,3,1))
plot(100*comparisons$observed,1:2,xlim=c(-25,40),ylim=c(.6,2.4),yaxt="n",ylab="",pch=1,
 xlab="Difference (percentage points)",main="Observed and interventional comparisons")
axis(2,at=1:2,labels=comparisons$case,las=1)
segments(100*comparisons$observed,1:2,100*comparisons$intervention,1:2,col="grey")
points(100*comparisons$intervention,1:2,pch=19,col="steelblue");abline(v=0,lty=2)
legend("topright",c("Observed","Intervention under known toy process"),pch=c(1,19),
 col=c("black","steelblue"),bty="n",cex=.8)
invisible(dev.off())
cat("END_41b\n")

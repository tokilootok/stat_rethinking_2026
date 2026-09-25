# A09 Problem 1B: base R, invented data. Run from repository root.

# Section 14: Reproduce the aggregate comparison
cat("BEGIN_SECTION_14\n")
# Invented data, not UCBadmit. Run the blocks in section order.
cells <- data.frame(category = c(1,1,2,2),
  department = c("Selective","Less selective","Selective","Less selective"),
  admitted = c(16,16,2,56), applied = c(80,20,20,80))
cells$rate <- with(cells, admitted/applied)
print(cells, row.names=FALSE)
totals <- aggregate(cbind(admitted,applied) ~ category, cells, sum)
totals$rate <- with(totals, admitted/applied)
print(totals, row.names=FALSE)
cat("Raw aggregate contrast (percentage points):",100*(totals$rate[1]-totals$rate[2]),"\n")
cat("END_SECTION_14\n")

# Section 15: Fit the stated logistic-normal models
cat("BEGIN_SECTION_15\n")
# Fit independent alpha ~ Normal(0,1) priors on logits, as in the text.
# Each parameter has a one-dimensional posterior, evaluated on a grid.
fit_grid <- function(y,n,step=.0025) {
  alpha <- seq(-12,12,by=step)
  log_mass <- dbinom(y,n,plogis(alpha),log=TRUE)+dnorm(alpha,0,1,log=TRUE)
  mass <- exp(log_mass-max(log_mass)); mass <- mass/sum(mass)
  stopifnot(sum(mass[abs(alpha)>11]) < 1e-8)
  list(p=plogis(alpha),mass=mass)
}
summarize <- function(x) c(mean=mean(x),lo90=unname(quantile(x,.05)),
                           hi90=unname(quantile(x,.95)))
set.seed(9015)
S <- 100000
cell_fits <- Map(fit_grid,cells$admitted,cells$applied)
p_draw <- vapply(cell_fits,function(f)
  sample(f$p,S,replace=TRUE,prob=f$mass),numeric(S))
# Columns retain the order of the four rows of cells.
within_draw <- cbind(Selective=p_draw[,1]-p_draw[,3],
                    Less_selective=p_draw[,2]-p_draw[,4])
cat("Within-department posterior contrasts, percentage points:\n")
print(round(100*t(apply(within_draw,2,summarize)),2))
aggregate_fits <- Map(fit_grid,totals$admitted,totals$applied)
p_aggregate <- vapply(aggregate_fits,function(f)
  sample(f$p,S,replace=TRUE,prob=f$mass),numeric(S))
aggregate_draw <- p_aggregate[,1]-p_aggregate[,2]
cat("Separate aggregate-model posterior contrast, percentage points:\n")
print(round(100*summarize(aggregate_draw),2))
# Verify posterior means against a grid twice as fine.
finer <- Map(function(y,n) fit_grid(y,n,.00125),cells$admitted,cells$applied)
mean_p <- function(f) sum(f$p*f$mass)
stopifnot(max(abs(vapply(cell_fits,mean_p,numeric(1))-
                 vapply(finer,mean_p,numeric(1)))) < 1e-6)
cat("END_SECTION_15\n")

# Section 16: Reconstruct the mixture in R
cat("BEGIN_SECTION_16\n")
cells$own_weight <- cells$applied/ave(cells$applied,cells$category,FUN=sum)
cells$contribution <- cells$rate*cells$own_weight
print(cells[c("category","department","rate","own_weight","contribution")],
      row.names=FALSE)
print(aggregate(contribution ~ category,cells,sum),row.names=FALSE)
stopifnot(isTRUE(all.equal(cells$rate[1:2]-cells$rate[3:4],c(.1,.1))))
cat("END_SECTION_16\n")

# Section 18: Standardize the same posterior draws
cat("BEGIN_SECTION_18\n")
departments <- cells$department[1:2]
pooled_n <- tapply(cells$applied,cells$department,sum)
w <- pooled_n[departments]/sum(pooled_n)
print(w)
cat("Raw standardized probabilities:",sum(w*cells$rate[1:2]),
    sum(w*cells$rate[3:4]),"\n")
standardized_draw <- drop(within_draw %*% as.numeric(w))
contrasts <- cbind(Aggregate=aggregate_draw,within_draw,Standardized=standardized_draw)
cat("Posterior contrasts, percentage points; central 90% intervals:\n")
print(round(100*t(apply(contrasts,2,summarize)),2))
cat("Pr(standardized contrast > 0):",round(mean(standardized_draw>0),3),"\n")
# Expanding labels verifies weighting; it adds no individual information.
d_id <- match(rep(cells$department,cells$applied),departments)
stopifnot(max(abs(rowMeans(within_draw[1:1000,d_id,drop=FALSE])-
                 standardized_draw[1:1000])) < 1e-12)
dir.create("courses/figures/A09_problem1b",recursive=TRUE,showWarnings=FALSE)
png("courses/figures/A09_problem1b/contrasts.png",width=1200,height=650,res=130)
z <- 100*t(apply(contrasts,2,summarize))
par(mar=c(5,9,3,1))
plot(z[,"mean"],1:4,xlim=range(z[,c("lo90","hi90")]),yaxt="n",ylab="",
     xlab="Category 1 minus category 2 (percentage points)",pch=19,
     main="Invented admissions data: posterior contrasts")
axis(2,at=1:4,labels=rownames(z),las=1)
segments(z[,"lo90"],1:4,z[,"hi90"],1:4,lwd=3,col="steelblue")
points(z[,"mean"],1:4,pch=19); abline(v=0,lty=2)
invisible(dev.off())
cat("END_SECTION_18\n")

# Section 19: Change the target weights explicitly
cat("BEGIN_SECTION_19\n")
targets <- rbind(Pooled=c(.5,.5),Mostly_selective=c(.8,.2),
                 Mostly_less_selective=c(.2,.8))
weighted_draw <- within_draw %*% t(targets)
print(round(100*t(apply(weighted_draw,2,summarize)),2))
cat("END_SECTION_19\n")

# Section 20: Check the calculation scale
cat("BEGIN_SECTION_20\n")
a1 <- qlogis(.20); a2 <- qlogis(.10)
print(round(c(log_odds_difference=a1-a2,odds_ratio=exp(a1-a2),
  probability_difference=plogis(a1)-plogis(a2),
  incorrect_difference=plogis(a1-a2)),4))
cat("END_SECTION_20\n")

# Section 21: Simulate future admission counts
cat("BEGIN_SECTION_21\n")
# New cohort: 80 category-1 applicants in the selective department.
set.seed(9021)
expected_admissions <- 80*p_draw[,1]
future_admissions <- rbinom(S,size=80,prob=p_draw[,1])
print(round(rbind(Expected_count=summarize(expected_admissions),
                  Future_observed_count=summarize(future_admissions)),2))
cat("END_SECTION_21\n")

# Section 22: Check replicated admissions in each cell
cat("BEGIN_SECTION_22\n")
set.seed(9022)
y_rep <- vapply(seq_len(nrow(cells)),function(j)
  rbinom(S,size=cells$applied[j],prob=p_draw[,j]),numeric(S))
predictive <- t(apply(y_rep,2,summarize))
print(cbind(cells[c("category","department","admitted")],
            round(predictive,2)),row.names=FALSE)
png("courses/figures/A09_problem1b/predictive_check.png",width=1200,height=650,res=130)
par(mar=c(6,5,3,1))
plot(1:4,predictive[,"mean"],ylim=c(0,max(predictive[,"hi90"])+3),
     xaxt="n",xlab="",ylab="Admissions",pch=1,
     main="Cell-level posterior predictive check (invented data)")
axis(1,at=1:4,labels=paste("Category",cells$category,
  ifelse(cells$department=="Selective","selective","less selective")),cex.axis=.8)
segments(1:4,predictive[,"lo90"],1:4,predictive[,"hi90"],col="steelblue",lwd=3)
points(1:4,predictive[,"mean"],pch=1); points(1:4,cells$admitted,pch=19)
legend("topleft",c("Observed","Predictive mean","Central 90% interval"),
  pch=c(19,1,NA),lty=c(NA,NA,1),col=c("black","black","steelblue"),bty="n")
invisible(dev.off())
cat("END_SECTION_22\n")

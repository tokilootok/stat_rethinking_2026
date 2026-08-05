# black cats new script

library(rethinking)
library(survival)

data(AustinCats)
d <- AustinCats

# prepare data for Stan models
dat <- list(
    N = nrow(d),
    days = d$days_to_event,
    adopted = ifelse( d$out_event=="Adoption" , 1 , 0 ),
    color = ifelse( d$color=="Black" , 1 , 2 ) )

# plot cats as survival lines
n <- 100
idx <- sample(1:dat$N,size=n)
ymax <- max(dat$days[idx])
plot(NULL,xlim=c(0,ymax),ylim=c(1,n),xlab="days observed",ylab="cat")
for ( i in 1:n ) { 
    j <- idx[i]
    cat_color <- ifelse( dat$color[j]==1 , "black" , "orange" )
    lines( c(0,dat$days[j]) , c(i,i) , lwd=4 , col=cat_color )
    if ( dat$adopted[j]==1 ) points( dat$days[j] , i , pch=16 , cex=1.5 , col=cat_color )
}

# generative code


cat_adopt <- function(day,prob) {
    if ( day > 1000 ) return(day)
    if ( runif(1) > prob ) {
        # keep waiting...
        day <- cat_adopt( day+1 , prob )
    }
    return( day ) # adopted
}

sim_cats1 <- function(n=10,p=c(0.1,0.2)) {
    color <- rep(NA,n)
    days <- rep(NA,n)
    for ( i in 1:n ) {
        color[i] <- sample(c(1,2),size=1,replace=TRUE)
        days[i] <- cat_adopt(1,p[color[i]])
    }
    return(list(N=n,days=days,color=color,adopted=rep(1,n)))
}

synth_cats <- sim_cats1(1e3)
## plot empirical K-M curves 
sfit <- survfit(Surv(days, adopted) ~ color, data = synth_cats)
plot( sfit, lty = 1 , lwd=3 , col=c("black","orange") , xlab="days" , ylab="proportion un-adopted" ) 

# model 1 - no censoring - just need geometric estimator

## prior predictive
n <- 12
sim_prior <- replicate(n,rbeta(2,1,10))
cols <- c( "black" , "orange" )
# rethinking::blank2(w=1.2)
plot( NULL , xlab="days" , ylab="proportion un-adopted" , xlim=c(0,50) , ylim=c(0,1) )
mtext("Prior predictive distribution")
for ( i in 1:n ) {
    days_rep <- sim_cats1(n=1e3,p=sim_prior[,i])
    xfit <- survfit(Surv(days, adopted) ~ color, data = days_rep)
    lines( xfit , lwd=2 , col=cols )
}


## test the first model code
p <- c(0.1,0.15)
# p <- c(0.15,0.1)
sim_dat <- sim_cats1(n=1000,p=p)

m1 <- ulam(
    alist(
        adopted|adopted==1 ~ bernoulli(p),
        p <- (1-P[color])^(days-1) * P[color],
        P[color] ~ beta(1,10)
    ), data=sim_dat , chains=4 , cores=4 )

postx <- extract.samples(m1)

plot(density( postx$P[,1] ) , lwd=3 , xlab="probability of adoption" , xlim=c(0.07,0.2) , main="" , ylim=c(0,100) )
k <- density( postx$P[,2] )
lines( k$x , k$y , lwd=3 , col="orange" )
abline(v=p[1],lwd=2)
abline(v=p[2],lwd=2,col="orange")

## real sample

m1 <- ulam(
    alist(
        adopted|adopted==1 ~ custom( log( (1-P[color[i]])^(days[i]-1) * P[color[i]] ) ),
        P[color] ~ beta(1,10)
    ), data=dat , chains=4 , cores=4 )

post1 <- extract.samples(m1)

# kaplan-meier posterior simulations
cols <- c( col.alpha("black") , col.alpha("orange") )
# rethinking::blank2(w=1.2)
plot( NULL , xlab="days" , ylab="proportion un-adopted" , xlim=c(0,50) , ylim=c(0,1) )
mtext("Posterior predictive distribution (1000 cats)")
n <- 12
for ( i in 1:n ) {
    days_rep <- sim_cats1(n=1e3,p=post1$P[i,])
    xfit <- survfit(Surv(days, adopted) ~ color, data = days_rep)
    lines( xfit , lwd=2 , col=cols )
}


## now add observation (censoring)

sim_cats2 <- function(n=10,p=c(0.1,0.2),cens=30) {
    color <- rep(NA,n)
    days <- rep(NA,n)
    for ( i in 1:n ) {
        color[i] <- sample(c(1,2),size=1,replace=TRUE)
        days[i] <- cat_adopt(1,p[color[i]])
    }
    adopted <- ifelse( days < cens , 1 , 0 )
    days <- ifelse( adopted==1 , days , cens )
    return(list(N=n,days=days,color=color,adopted=adopted))
}

## test censoring model
p <- c(0.1,0.15)
sim_dat <- sim_cats2(n=1e3,p=p,cens=20)

plot( sim_dat$days , xlab="cat" , ylab="days (simulated)" , col=ifelse(sim_dat$color==1,1,"orange") )

m2 <- ulam(
    alist(
        # adopted
        adopted|adopted==1 ~ custom( log( (1-P[color[i]])^(days[i]-1) * P[color[i]] ) ),
        # censored
        adopted|adopted==0 ~ custom( log( (1-P[color[i]])^days[i] ) ),
        P[color] ~ beta(1,10)
    ), data=sim_dat , chains=4 , cores=4 , sample=TRUE )

# compare to model that ignores censoring
m1 <- ulam(
    alist(
        adopted|adopted==1 ~ custom( log( (1-P[color[i]])^(days[i]-1) * P[color[i]] ) ),
        P[color] ~ beta(1,10)
    ), data=sim_dat , chains=4 , cores=4 )

post1 <- extract.samples(m1)
post2 <- extract.samples(m2)

plot(density( post2$P[,1] ) , lwd=3 , xlab="probability of adoption" , xlim=c(0.07,0.2) , main="" , ylim=c(0,100) )
k <- density( post2$P[,2] )
lines( k$x , k$y , lwd=3 , col="orange" )
abline(v=p[1],lwd=2)
abline(v=p[2],lwd=2,col="orange")

k <- density( post1$P[,1] )
lines( k$x , k$y , lwd=3 , col="black" , lty=2 )
k <- density( post1$P[,2] )
lines( k$x , k$y , lwd=3 , col="orange" , lty=2 )

## now real sample

m2 <- ulam(
    alist(
        # adopted
        adopted|adopted==1 ~ custom( log( (1-P[color[i]])^(days[i]-1) * P[color[i]] ) ),
        # censored
        adopted|adopted==0 ~ custom( log( (1-P[color[i]])^days[i] ) ),
        P[color] ~ beta(1,10)
    ), data=dat , chains=4 , cores=4 , sample=TRUE )

# will compare to inferences from model that ignores censoring
m1 <- ulam(
    alist(
        adopted|adopted==1 ~ custom( log( (1-P[color[i]])^(days[i]-1) * P[color[i]] ) ),
        P[color] ~ beta(1,10)
    ), data=dat , chains=4 , cores=4 )

post1 <- extract.samples(m1)
post2 <- extract.samples(m2)

# kaplan-meier posterior simulations
cols <- c( col.alpha("black") , col.alpha("orange") )
# rethinking::blank2(w=1.2)
plot( NULL , xlab="days" , ylab="proportion un-adopted" , xlim=c(0,50) , ylim=c(0,1) )
mtext("Posterior predictive distribution (1000 cats)")
n <- 12
for ( i in 1:n ) {
    days_rep <- sim_cats1(n=1e3,p=post1$P[i,])
    xfit <- survfit(Surv(days, adopted) ~ color, data = days_rep)
    lines( xfit , lwd=2 , col=cols )
}

for ( i in 1:n ) {
    days_rep <- sim_cats1(n=1e3,p=post2$P[i,])
    xfit <- survfit(Surv(days, adopted) ~ color, data = days_rep)
    lines( xfit , lwd=2 , col=cols )
}


## model with a parameter for each censored observation

cat_code3 <- "
// imputation version
data{
    int N;
    array[N] int adopted; // 1/0 indicator
    vector[N] days;       // days until event
    array[N] int color;   // 1=black, 2=other
}
parameters{
    vector<lower=0,upper=1>[2] P;
    vector<lower=days>[N] days_imputed;
}
model{
    P ~ beta(1,10);
    for ( i in 1:N ) {
        real PC = P[color[i]];
        if ( adopted[i]==1 ) {
            target += log( (1-PC)^(days[i]-1) * PC );
            // just copy observed value into parameter vector
            days_imputed[i] ~ normal(days[i],0.01);
        } else {
            // same likelihood is now *prior* for missing value
            target += log( (1-PC)^(days_imputed[i]-1) * PC );
        }
    }
}
"

m3 <- cstan(model_code=cat_code3,data=dat,chains=4,cores=4)

precis(m3,2,pars="P")

plot(precis(m3,2,pars="days_imputed"))

m3u <- ulam(
    alist(
        adopted|adopted==1 ~ custom( log( (1-P[color[i]])^(days[i]-1) * P[color[i]] ) ),
        adopted|adopted==0 ~ custom( log( (1-P[color[i]])^(days_impute[i]-1) * P[color[i]] ) ),
        vector[N]:days_impute ~ uniform(days,2000),
        P[color] ~ beta(1,10)
    ), data=dat , chains=4 , cores=4 , sample=TRUE )


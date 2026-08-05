library(rethinking)

set.seed(12)
N <- 2000 # number of applicants
# even gender distribution
G <- sample( 1:2 , size=N , replace=TRUE )
# sample ability, high (1) to average (0)
u <- rbern(N,0.1)
# gender 1 tends to apply to department 1, 2 to 2
# and G=1 with greater ability tend to apply to 2 as well
D <- rbern( N , ifelse( G==1 , u*1 , 0.75 ) ) + 1
# matrix of acceptance rates [dept,gender]
p_u0 <- matrix( c(0.1,0.1,0.1,0.3) , nrow=2 )
p_u1 <- matrix( c(0.3,0.3,0.5,0.5) , nrow=2 )
p_u <- list( p_u0 , p_u1 )
# simulate acceptance
p <- sapply( 1:N , function(i) p_u[[1+u[i]]][D[i],G[i]] )
A <- rbern( N , p )

dat_sim <- list( A=A , D=D , G=G )

# total effect gender
m1 <- ulam(
    alist(
        A ~ bernoulli(p),
        logit(p) <- a[G],
        a[G] ~ normal(0,1)
    ), data=dat_sim , chains=4 , cores=4 )

# direct effects - now confounded!
m2 <- ulam(
    alist(
        A ~ bernoulli(p),
        logit(p) <- a[G,D],
        matrix[G,D]:a ~ normal(0,1)
    ), data=dat_sim , chains=4 , cores=4 )

post2 <- extract.samples(m2)
dens(inv_logit(post2$a[,1,1]),lwd=3,col=4,xlim=c(0,0.5),xlab="probability admission")
dens(inv_logit(post2$a[,2,1]),lwd=3,col=4,lty=2,add=TRUE)
dens(inv_logit(post2$a[,1,2]),lwd=3,col=2,add=TRUE)
dens(inv_logit(post2$a[,2,2]),lwd=3,col=2,lty=2,add=TRUE)

# contrasts

post2 <- extract.samples(m2)

post2$fm_contrast_D1 <- 
    post2$a[,1,1] - post2$a[,2,1]

post2$fm_contrast_D2 <- 
    post2$a[,1,2] - post2$a[,2,2]

dens(post2$fm_contrast_D1)
dens(post2$fm_contrast_D2,add=TRUE,col=2)

# sensitivity

datl <- dat_sim
datl$D2 <- ifelse(datl$D==2,1,0)
datl$N <- length(datl$D)
datl$b <- c(1,1)
datl$g <- c(1,0)

mGDu <- ulam(
    alist( 
        # A model
        A ~ bernoulli(p),
        logit(p) <- a[G,D] + b[G]*u[i],
        matrix[G,D]:a ~ normal(0,1),

        # D model
        D2 ~ bernoulli(q),
        logit(q) <- delta[G] + g[G]*u[i],
        delta[G] ~ normal(0,1),

        # declare unobserved u
        vector[N]:u ~ normal(0,1)
    ), data=datl , chains=4 , cores=4 )

post2 <- extract.samples(mGDu)

# now treat b and g vectors as priors instead of fixed values

datl$b <- NULL
datl$g <- NULL
mGDu2 <- ulam(
    alist( 
        # A model
        A ~ bernoulli(p),
        logit(p) <- a[G,D] + b[G]*u[i],
        matrix[G,D]:a ~ normal(0,1),
        vector[2]:b ~ uniform(0,1),

        # D model
        D2 ~ bernoulli(q),
        logit(q) <- delta[G] + g[G]*u[i],
        delta[G] ~ normal(0,1),
        vector[2]:g ~ uniform(0,1),

        # declare unobserved u
        vector[N]:u ~ normal(0,1)
    ), data=datl , chains=4 , cores=4 )

post2 <- extract.samples(mGDu2)

dens(inv_logit(post2$a[,1,1]),lwd=3,col=4,xlim=c(0,0.5),xlab="probability admission")
dens(inv_logit(post2$a[,2,1]),lwd=3,col=4,lty=2,add=TRUE)
dens(inv_logit(post2$a[,1,2]),lwd=3,col=2,add=TRUE)
dens(inv_logit(post2$a[,2,2]),lwd=3,col=2,lty=2,add=TRUE)

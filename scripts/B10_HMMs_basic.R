
library(rethinking)

# plain vanilla Stan HMM

basic_hmm_code <- "
data {
  int N;            // Number of days
  int K;            // Number of states: 1=alive, 2=dead
  array[N] int Y;   // 0/1 President seen on day i
}
parameters {
  // probability of death per day
  real<lower=0,upper=1> phi;

  // probability seen conditional on state
  vector<lower=0,upper=1>[K] theta;
}
transformed parameters {
  // Build transition matrix
  matrix[K, K] gamma;
  gamma[1,1] = 1-phi; // alive->alive
  gamma[1,2] = phi;   // alive->dead
  gamma[2,1] = 0;     // dead->alive
  gamma[2,2] = 1;     // dead->dead

  // Initial state
  vector[K] rho = [1.0,0.0]'; // Alive at start

  // Compute the log likelihoods in each possible state
  matrix[K, N] log_omega;
  for (n in 1:N) {
    for ( k in 1:K )
      log_omega[k, n] = bernoulli_lpmf( Y[n] | theta[k] );
  }
}
model {
  // prior for emissions for each state
  theta[1] ~ beta(1,4);
  theta[2] ~ beta(1,100); // not zero exactly
  
  // mortality prior
  phi ~ beta(1,50); // mean prob 1/(50+1) ~= 0.02

  // Increment target by log p(Y | phi, gamma, rho)
  target += hmm_marginal(log_omega, gamma, rho);
}
generated quantities {
  array[N] int latent_states = hmm_latent_rng(log_omega, gamma, rho);
  matrix[K, N] hidden_probs = hmm_hidden_state_prob(log_omega, gamma, rho);
}
"

# example data
library(rethinking)

# number of days
N <- 60
# simulate death
X <- rep(1,N)
( dead_after <- rgeom(1,0.05) )
dead_after <- 50 # dies on day 50
if (dead_after < N) X[dead_after:N] <- 2
# simulate public appearances
Y <- rep(0,N)
for ( i in 1:N ) {
    if ( X[i]==1 ) # alive, appear 50% of days
        Y[i] <- rbern(1,0.5) 
}

# plot data
blank2(w=1.7)

plot(X,xlab="Day",ylab="State",lwd=2,col=2)
abline(v=dead_after,col=2,lty=2)

plot(Y,xlab="Day",ylab="Appeared in public",lwd=2,col=2)
abline(v=dead_after,col=2,lty=2)

# analyze
dat <- list(
    Y = Y,
    N = N,
    K = 2
)

m0 <- cstan(model_code=basic_hmm_code,data=dat,chains=4,core=4)

precis(m0,2)

# blank2(w=1.7)

post <- extract.samples(m0)

pd <- post$hidden_probs[,2,]
plot( 1:N , apply(pd,2,mean) , col=2 , lwd=2 , xlab="Day" , ylab="Probability Dead" , ylim=c(0,1) )
pi <- apply(pd,2,PI)
for ( i in 1:N ) lines(c(i,i),pi[,i],col=2)

apply(post$gamma,2:3,mean)

# now revise example for many presidents - actually a capture-recapture model
# we want to estimate mortality from annual surveys

cr_hmm_code <- "
data {
  int N;            // Number of years
  int M;            // Number of individuals
  int K;            // Number of states: 1=alive, 2=dead
  array[M,N] int Y;   // 0/1 individual n seen on year m
}
parameters {
  // probability of death per year
  real<lower=0,upper=1> phi;

  // probability seen conditional on state
  vector<lower=0,upper=1>[K] theta;
}
transformed parameters {
  // Build transition matrix
  matrix[K, K] gamma;
  gamma[1,1] = 1-phi; // alive->alive
  gamma[1,2] = phi;   // alive->dead
  gamma[2,1] = 0;     // dead->alive
  gamma[2,2] = 1;     // dead->dead

  // Initial state
  vector[K] rho = [1.0,0.0]'; // Alive at start (first capture)

  // Compute the log likelihoods in each possible state
  array[M] matrix[K, N] log_omega;
  for ( m in 1:M )      // individuals
    for (n in 1:N) {    // observations
      for ( k in 1:K )  // states
          log_omega[m, k, n] = bernoulli_lpmf( Y[m,n] | theta[k] );
    }
}
model {
  // prior for emissions for each state
  theta[1] ~ beta(1,4);
  theta[2] ~ beta(1,100); // not zero exactly
  
  // mortality prior
  phi ~ beta(1,10); // mean prob 1/(10+1) ~= 0.09

  // Increment target by log p(Y | phi, gamma, rho)
  for ( m in 1:M )
    target += hmm_marginal(log_omega[m], gamma, rho);
}
generated quantities {
  //array[N] int latent_states = hmm_latent_rng(log_omega, gamma, rho);
  //matrix[K, N] hidden_probs = hmm_hidden_state_prob(log_omega, gamma, rho);
}
"

# simulate some records

# number of years, individuals
N <- 20
M <- 100
Y <- matrix(0,nrow=M,ncol=N)
X <- matrix(1,nrow=M,ncol=N)
# prob death per year
phi <- 0.05
# simulate deaths
for ( m in 1:M ) {
  dead_after <- rgeom(1,phi) + 2
  if (dead_after < N) X[m,dead_after:N] <- 2
}
# simulate captures
for ( m in 1:M )
  for ( i in 1:N ) {
      if ( X[m,i]==1 ) # alive, appear 50% of surveys
          Y[m,i] <- rbern(1,0.5) 
  }

dat2 <- list(
  N=N,
  M=M,
  Y=Y,
  K=2
)

m1 <- cstan(model_code=cr_hmm_code,data=dat2,chains=4,core=4)

precis(m1,2)

post <- extract.samples(m1)
# blank2()
dens(post$phi,col=2,xlab="phi (annual mortality)",xlim=c(0,0.2),lwd=2)
curve(dbeta(x,1,10),add=TRUE,lty=2,lwd=2)
abline(v=phi,lty=2,col=1)

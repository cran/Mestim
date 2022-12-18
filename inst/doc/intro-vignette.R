## ---- message=FALSE, warning=FALSE--------------------------------------------
gen_lr_dat <- function(n, seed=123)
{
set.seed(seed)
X_1 <- rnorm(n, sd = 1); X_2 <- rnorm(n, sd = 3) # generate x_1 and x_2
true_betas <- c(4,5) # generate true parameters
X <- model.matrix(~-1+X_1+X_2) # build the design matrix
Y <- rbinom(n, 1, 1/(1 + exp(-X %*% true_betas)) ) # generate Y from X and true_betas
dat  <-  data.frame(X_1=X_1, X_2=X_2, Y=Y) # build a simulated dataset
return(dat)
}

## ---- message=FALSE, warning=FALSE--------------------------------------------
dat <- gen_lr_dat(5000)
head(dat)

## ---- message=FALSE, warning=FALSE, results='hide'----------------------------
mod <- glm(Y~-1 + X_1 + X_2, data=dat, family = "binomial")
thetas_hat <- list(thetas_1=coef(mod)[1], thetas_2=coef(mod)[2])

## ---- message=FALSE, warning=FALSE, results='hide'----------------------------
psi_1 <- expression( ((1/(1+exp(-(thetas_1 * X_1 + thetas_2 * X_2)))) - Y) * X_1 )
psi_2 <- expression( ((1/(1+exp(-(thetas_1 * X_1 + thetas_2 * X_2)))) - Y) * X_2 )
psi <- list(psi_1, psi_2)

## ---- message=FALSE, warning=FALSE, results='hide'----------------------------
library(Mestim)
res <- get_vcov(data=dat, thetas=thetas_hat, M=psi)

## ---- message=FALSE, warning=FALSE--------------------------------------------
res$vcov

## ---- message=FALSE, warning=FALSE--------------------------------------------
vcov(mod)

## ---- message=FALSE, warning=FALSE--------------------------------------------
res$se

## ---- message=FALSE, warning=FALSE--------------------------------------------
gen_obs_dat <- function(n, seed=123)
{
set.seed(seed)
X <- rnorm(n) # generate X
A <- rbinom(n, 1, 1/(1 + exp(- 2 * X)) )  # generate treatment allocation A
X_mat <- model.matrix(~ -1 + X + A + A:X) # build the design matrix
true_gammas <- c(4,3,2) 
epsi <- rnorm(n,0,20) # generate gaussian noise 
Y <- (X_mat %*% true_gammas) + epsi # generate observed outcomes 
dat  <-  data.frame(X=X, A=A, Y=Y) # build a simulated dataset
return(dat)
}

## ---- message=FALSE, warning=FALSE--------------------------------------------
dat <- gen_obs_dat(5000)
head(dat)

## ---- message=FALSE, warning=FALSE, results='hide'----------------------------
m <- lm(Y~ -1 + X + A + A:X, data = dat)
gamma_1_hat <- coef(m)[1]
gamma_2_hat <- coef(m)[2]
gamma_3_hat <- coef(m)[3]

## ---- message=FALSE, warning=FALSE, results='hide'----------------------------
delta_hat <- mean(gamma_2_hat + gamma_3_hat*dat$X)

## ---- message=FALSE, warning=FALSE, results='hide'----------------------------
psi_1 <- expression( (Y - gamma_1*X - gamma_2*A - gamma_3*A*X) * X )
psi_2 <- expression( (Y - gamma_1*X - gamma_2*A - gamma_3*A*X) * A )
psi_3 <- expression( (Y - gamma_1*X - gamma_2*A - gamma_3*A*X) * A*X )
psi_4 <- expression( gamma_2 + gamma_3 * X - delta )
psi <- list(psi_1, psi_2, psi_3, psi_4)

## ---- message=FALSE, warning=FALSE, results='hide'----------------------------
thetas_hat <- list(gamma_1=gamma_1_hat,
                  gamma_2=gamma_2_hat,
                  gamma_3=gamma_3_hat,
                  delta=delta_hat)

## ---- message=FALSE, warning=FALSE--------------------------------------------
res <- get_vcov(data=dat, thetas=thetas_hat, M=psi)
res$vcov

## ---- message=FALSE, warning=FALSE--------------------------------------------
boot_fun <- function(d, i=1:nrow(d)) {
  z<-d[i,]
  mod <- lm(Y~ -1 + X + A + A:X, data = z)
  gamma_1_hat <- coef(mod)[1]
  gamma_2_hat <- coef(mod)[2]
  gamma_3_hat <- coef(mod)[3]
  delta_hat <- mean(gamma_2_hat*1 + gamma_3_hat*1*z$X)
  return( c(gamma_1_hat, gamma_2_hat, gamma_3_hat, delta_hat) )
}
boot_start_time <- Sys.time()
res_boot <- boot::boot(dat, boot_fun, R=999)
boot_end_time <- Sys.time()
paste("Bootstrapping took", round(as.numeric(boot_end_time - boot_start_time), 2), "seconds.")

## ---- message=FALSE, warning=FALSE, echo=FALSE--------------------------------
cov(res_boot$t)
Mestim_start_time <- Sys.time()
res <- get_vcov(data=dat, thetas=thetas_hat, M=psi)
Mestim_end_time <- Sys.time()
#paste("Mestim took", round(as.numeric(Mestim_end_time - Mestim_start_time), 2), "seconds.")

## ---- message=FALSE, warning=FALSE--------------------------------------------
gen_dtr_dat <- function(n, seed=456)
{
set.seed(seed)
expit <- function(x) 1/(1+exp(-x))

X_1 <- rnorm(n, sd=.1)
S_1 <- exp(rnorm(n, mean = X_1, sd = .1))
A_1 <- rbinom(n, size = 1, prob = expit(-.1+log(S_1)))

X_2 <- (X_1>0) * rnorm(n, mean = 1.1*X_1 - .5 * A_1, sd=.05) + (X_1<0) * X_1
S_2 <- exp(rnorm(n, mean = X_2, sd = .1))
A_2 <- rbinom(n, size = 1, prob = expit(.1+log(S_2)+3*A_1))

X_3 <- (X_2>0) * rnorm(n, mean = 1.1*X_2 - .5 * A_2, sd=.05) + (X_2<0) * X_2
Y <- exp(rnorm(n, mean = X_3 + .1*(A_1 + A_2), sd = .1)) #0.1 penalty for treating

dat <- data.frame(S_1=S_1, A_1=A_1, S_2=S_2, A_2=A_2, Y)  
return(dat)
}

## ---- message=FALSE, warning=FALSE--------------------------------------------
dat <- gen_dtr_dat(5000)
head(dat)

## ---- message=FALSE, warning=FALSE, results='hide'----------------------------
e_1 <- glm(A_1~I(log(S_1)), data=dat, family = "binomial")
delta_1_hat <- coef(e_1)[1]
delta_2_hat <- coef(e_1)[2]

e_2 <- glm(A_2~I(log(S_2)) + A_1 , data=dat, family = "binomial")
phi_1_hat <- coef(e_2)[1]
phi_2_hat <- coef(e_2)[2]
phi_3_hat <- coef(e_2)[3]

## ---- message=FALSE, warning=FALSE, results='hide'----------------------------
dat$log_S_1 <- log(dat$S_1) ; dat$log_S_2 <- log(dat$S_2) # For ease of programming

psi_1 <- expression( ((1/(1+exp(-(delta_1 + delta_2 * log_S_1)))) - A_1) * 1 )
psi_2 <- expression( ((1/(1+exp(-(delta_1 + delta_2 * log_S_1)))) - A_1) * log_S_1)

psi_3 <- expression( ((1/(1+exp(-(phi_1 + phi_2 * log_S_2 + phi_3 * A_1)))) - A_2) * 1 )
psi_4 <- expression( ((1/(1+exp(-(phi_1 + phi_2 * log_S_2 + phi_3 * A_1)))) - A_2) * log_S_2)
psi_5 <- expression( ((1/(1+exp(-(phi_1 + phi_2 * log_S_2 + phi_3 * A_1)))) - A_2) * A_1)

## ---- message=FALSE, warning=FALSE, results='hide'----------------------------
# The regime we are interested in
dat$d_1 <- dat$S_1>1
dat$d_2 <- dat$S_2>1

# For ease of programming
dat$C_d <- with(dat, d_1==A_1 & d_2==A_2)

# Store the last element of psi
psi_6 <- expression( Y * C_d *
          (  (1+exp(-(delta_1 + delta_2 * log_S_1)))^d_1 * # numerator
             (1+exp(-(phi_1 + phi_2 * log_S_2 + phi_3 * A_1)))^d_2  ) 
          
      /   (  (1-(1+exp(-(delta_1 + delta_2 * log_S_1)))^(-1))^(1-d_1) * # denominator
             (1-(1+exp(-(phi_1 + phi_2 * log_S_2 + phi_3 * A_1)))^(-1))^(1-d_2)  ) 
- V
)

## ---- message=FALSE, warning=FALSE, results='hide'----------------------------
psi <- list(psi_1, psi_2, psi_3, psi_4, psi_5, psi_6)

## ---- message=FALSE, warning=FALSE, results='hide'----------------------------
# Just delete "- V" from in the previous expression
# add _hat as appropriate
ipw_estimator <- expression( Y * C_d *
          (  (1+exp(-(delta_1_hat + delta_2_hat * log_S_1)))^d_1 * # numerator
             (1+exp(-(phi_1_hat + phi_2_hat * log_S_2 + phi_3_hat * A_1)))^d_2  ) 
          
      /   (  (1-(1+exp(-(delta_1_hat + delta_2_hat * log_S_1)))^(-1))^(1-d_1) * # denominator
             (1-(1+exp(-(phi_1_hat + phi_2_hat * log_S_2 + phi_3_hat * A_1)))^(-1))^(1-d_2)  )  
)

# Compute individual contribution and take the average
V_hat <- with(dat, mean(eval(ipw_estimator))) # Other ways to compute this quantity are OK too.

thetas_hat <- list(delta_1=delta_1_hat,
                   delta_2=delta_2_hat,
                   phi_1=phi_1_hat,
                   phi_2=phi_2_hat,
                   phi_3=phi_3_hat,
                   V=V_hat)

## ---- message=FALSE, warning=FALSE--------------------------------------------
res <- get_vcov(data=dat, thetas=thetas_hat, M=psi)
res$se

## ---- message=FALSE, warning=FALSE--------------------------------------------
boot_fun <- function(d, i=1:nrow(d)) {
  z<-d[i,]
  e_1 <- glm(A_1~I(log(S_1)), data=z, family = "binomial")
  e_2 <- glm(A_2~I(log(S_2)) + A_1 , data=z, family = "binomial")

  delta_1_hat <- coef(e_1)[1]
  delta_2_hat <- coef(e_1)[2]
  phi_1_hat <- coef(e_2)[1]
  phi_2_hat <- coef(e_2)[2]
  phi_3_hat <- coef(e_2)[3]

  ipw_estimator <- expression( z$Y * z$C_d *
            (  (1+exp(-(delta_1_hat + delta_2_hat * z$log_S_1)))^z$d_1 * # numerator
               (1+exp(-(phi_1_hat + phi_2_hat * z$log_S_2 + phi_3_hat * z$A_1)))^z$d_2  ) 
            
        /   (  (1-(1+exp(-(delta_1_hat + delta_2_hat * z$log_S_1)))^(-1))^(1-z$d_1) * # denominator
               (1-(1+exp(-(phi_1_hat + phi_2_hat * z$log_S_2 + phi_3_hat * z$A_1)))^(-1))^(1-z$d_2)  )  
  )

  V_hat <- mean(eval(ipw_estimator))
  return( c(delta_1_hat, delta_2_hat, phi_1_hat, phi_2_hat, phi_3_hat, V_hat) )
}
boot_start_time <- Sys.time()
res_boot <- boot::boot(dat, boot_fun, R=999)
boot_end_time <- Sys.time()
paste("Bootstrapping took", round(as.numeric(boot_end_time - boot_start_time), 2), "seconds.")

## ---- message=FALSE, warning=FALSE, echo=FALSE--------------------------------
res_boot
Mestim_start_time <- Sys.time()
res <- get_vcov(data=dat, thetas=thetas_hat, M=psi)
Mestim_end_time <- Sys.time()
#paste("Mestim took", round(as.numeric(Mestim_end_time - Mestim_start_time), 2), "seconds.")


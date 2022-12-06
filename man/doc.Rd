\name{Mestim}
\alias{Mestim}
\alias{get_vcov}
\title{Parameters Variance-Covariance Matrix From M-estimation}
\description{
Provides a flexible framework for estimating the variance-covariance matrix of a set of estimated parameters.
Estimation relies on providing unbiased estimating functions to compute the empirical sandwich variance.
(i.e., M-estimation in the vein of Stefanski & Boos (2002) <doi:10.1198/000313002753631330>).
}
\usage{
get_vcov(data, betas, M)
}
\arguments{
  \item{data}{a dataframe with proper variable (i.e., column) names.}
  \item{betas}{a list of the (properly named) estimated parameters.}
  \item{M}{a list of expressions detailing the unbiased estimating functions with the same ordering as betas.
  The variables and varameters names in the expressions need be consistent with those of data and betas}
}
\value{
A list with elements \code{vcov}, \code{se}, and \code{jacob}.
\item{vcov}{pxp matrix of the estimated asymptotic variance-covariance matrix of the estimated parameters betas.}
\item{se}{length p vector of the estimated asymptotic standard error for the estimated parameters betas.}
\item{jacob}{a list of lists containing expressions for computing the jacobian matrix.}
}
\author{
François Grolleau <francois.grolleau@aphp.fr>
}
\references{
Stefanski, LA. and Boos DD. (2002)
\emph{The Calculus of M-Estimation, The American Statistician},
\doi{10.1198/000313002753631330}.
}
\examples{
####
## Simulate data
####
n <- 10000 # number of simulated iid observations
x_1 <- rnorm(n); x_2 <- rnorm(n) # generate x_1 and x_2
true_betas <- c(2,3) # generate true parameters
X <- model.matrix(~-1+x_1+x_2) # build the design matrix
y <- rbinom(n, 1, 1/(1 + exp(-X \%*\% true_betas)) ) # generate Y from X and true_betas
dat  <-  data.frame(x_1=x_1, x_2=x_2, y=y) # build a simulated dataset

####
## Fit a LR model (estimated parameters solve unbiased estimating equations)
####

mod <- glm(y~-1 + x_1 + x_2, data=dat, family = "binomial")

####
## Get variance covariance matrix for all parameters solving unbiased estimating equations
####

# Put estimated parameters in a list
betas_hat <- list(beta_1=coef(mod)[1], beta_2=coef(mod)[2])

# Build a list of unbiased estimating functions
# NB: parameters' names must be consistent with the previous list
M_1 <- expression( ((1/(1+exp( -( beta_1 * x_1 + beta_2 * x_2 ) ))) - y ) * x_1 )
M_2 <- expression( ((1/(1+exp( -( beta_1 * x_1 + beta_2 * x_2 ) ))) - y ) * x_2 )
est_functions <- list(M_1, M_2)

## Pass arguments and run get_vcov
res <- get_vcov(data=dat, betas=betas_hat, M=est_functions)

# Estimted variance covariance matrix is similar to that obtain from glm
res$vcov
vcov(mod)

# So are the standard errors for the estimated parameters
res$se
summary(mod)$coefficients[,2]

}

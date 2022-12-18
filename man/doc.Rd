\name{Mestim}
\alias{Mestim}
\alias{get_vcov}
\title{Parameters Variance-Covariance Matrix From M-estimation}
\description{
Provides a flexible framework for estimating the variance-covariance matrix of a multidimensional parameter.
Estimation relies on providing unbiased estimating functions to compute the empirical sandwich variance.
(i.e., M-estimation in the vein of Tsiatis et al. (2019) <doi:10.1201/9780429192692>).
}
\usage{
get_vcov(data, thetas, M)
}
\arguments{
  \item{data}{a dataframe with proper variable (i.e., column) names.}
  \item{thetas}{a list of the (properly named) estimated parameters.}
  \item{M}{a list of expressions detailing the unbiased estimating functions with the same ordering as \code{thetas}.
  The variables and parameters names in the expressions need be consistent with those of \code{data} and \code{thetas}.}
}
\value{
A list with elements \code{vcov}, \code{se}, and \code{jacob}.
\item{vcov}{pxp matrix of the estimated asymptotic variance-covariance matrix of the estimated parameters in \code{thetas}.}
\item{se}{length p vector of the estimated asymptotic standard error for the estimated parameters in \code{thetas}.}
\item{jacob}{a list of lists containing expressions for computing the jacobian matrix.}
}
\author{
François Grolleau <francois.grolleau@aphp.fr>
}
\references{
Stefanski, LA. and Boos DD. (2002)
\emph{The Calculus of M-Estimation, The American Statistician},
\doi{10.1198/000313002753631330}.\cr

Tsiatis, A. A., Davidian, M., Holloway, S. T. and Laber, E. B (2019)
\emph{Dynamic Treatment Regimes: Statistical Methods for Precision Medicine, CRC Press},
\doi{10.1201/9780429192692}.
}
\examples{
####
## Simulate data
####
set.seed(123)
n <- 10000 # number of simulated iid observations
x_1 <- rnorm(n); x_2 <- rnorm(n) # generate x_1 and x_2
true_thetas <- c(2,3) # generate true parameters
X <- model.matrix(~-1+x_1+x_2) # build the design matrix
y <- rbinom(n, 1, 1/(1 + exp(-X \%*\% true_thetas)) ) # generate Y from X and true_thetas
dat  <-  data.frame(x_1=x_1, x_2=x_2, y=y) # build a simulated dataset

####
## Fit a LR model (estimated parameters solve unbiased estimating equations)
####

mod <- glm(y~-1 + x_1 + x_2, data=dat, family = "binomial")

####
## Get variance covariance matrix for all parameters solving unbiased estimating equations
####

# Put estimated parameters in a list
thetas_hat <- list(theta_1=coef(mod)[1], theta_2=coef(mod)[2])

# Build a list of unbiased estimating functions
# NB: parameters' names must be consistent with the previous list
psi_1 <- expression( ((1/(1+exp( -( theta_1 * x_1 + theta_2 * x_2 ) ))) - y ) * x_1 )
psi_2 <- expression( ((1/(1+exp( -( theta_1 * x_1 + theta_2 * x_2 ) ))) - y ) * x_2 )
est_functions <- list(psi_1, psi_2)

## Pass arguments and run get_vcov
res <- get_vcov(data=dat, thetas=thetas_hat, M=est_functions)

# Estimted variance covariance matrix is similar to that obtain from glm
res$vcov
vcov(mod)

# So are the standard errors for the estimated parameters
res$se
summary(mod)$coefficients[,2]

}

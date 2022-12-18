get_vcov <- function (data, thetas, M)
{
  n <- nrow(data)
  list2env(data, environment())
  list2env(thetas, environment())
  thetas_names <- names(thetas)
  multi_deriv <- function(M_i) {
    lapply(thetas_names, function(beta) stats::D(M_i, beta))
  }
  jacob <- lapply(M, multi_deriv)
  mid <- sapply(M, function(it) eval(it))
  mid <- as.list(data.frame(t(mid)))
  mid <- lapply(mid, function(x) x %o% x)
  mid <- Reduce("+", mid)/n

  correct_len <- function(x) {
    if(length(x)==1){
    res <- rep(x, n)
    } else { res <- x}
    return(res)
    }

  comp_jacob <- function(x) {
    sapply(x, function(it) correct_len(eval(it)) )
  }
  lat <- sapply(jacob, comp_jacob, simplify = "array")
  lat <- lapply(seq(n), function(x) t(lat[x, , ]) )
  lat <- Reduce("+", lat)/n
  lat_inv <- solve(lat)
  Sigma <- lat_inv %*% mid %*% t(lat_inv)
  vcov <- Sigma/n
  colnames(vcov) <- thetas_names
  rownames(vcov) <- thetas_names
  se <- sqrt(diag(vcov))
  names(se) <- thetas_names
  res <- list(vcov = vcov, se = se, jacob = jacob)
  return(res)
}

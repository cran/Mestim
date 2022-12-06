get_vcov <- function(data, betas, M){

  n <- nrow(data)
  list2env(data, environment())
  list2env(betas, environment())
  betas_names <- names(betas)

  multi_deriv <- function(M_i){
    lapply(betas_names, function(beta) stats::D(M_i, beta))
  }

  Jacob <- lapply(M, multi_deriv)

  mid <- sapply(M, function(it) eval(it))
  mid <- as.list(data.frame(t(mid)))
  mid <- lapply(mid, function(x) x %o% x)
  mid <- Reduce('+', mid)/n

  comp_jacob <- function(x){
    sapply(x, function(it) eval(it))
  }

  lat <- sapply(Jacob, comp_jacob, simplify = "array")
  lat <- lapply(seq(dim(lat)[1]), function(x) lat[x, , ])
  lat <- Reduce('+', lat)/n

  lat_inv <- solve(lat)

  Sigma <- lat_inv %*% mid %*% t(lat_inv)
  vcov <- Sigma/n
  colnames(vcov) <- betas_names
  rownames(vcov) <- betas_names
  se <- sqrt(diag(vcov))
  names(se) <- betas_names

  res <- list(vcov=vcov, se=se, jacob=Jacob)

  return (  res  )

}

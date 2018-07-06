
### Bernoulli Mixture Model in Pure R
## Coded for understanding purposes
## Need to translate to C

make_test_data <- function(n, prototypes, pis) {

  samples <- list()
  for (i in seq.int(n)) {

    pi <- sample(seq_along(prototypes), size = 1, prob = pis/sum(pis))

    samples[[i]] <- sapply(prototypes[[pi]], function(p) rbinom(1, 1, p))

  }

  do.call(rbind, samples)
}

### code bernoulli mixture model in R

## Probability that xn was generated given a prototype
p_xn_i <- function(xn, proto) {
  lik <- 1
  for (i in seq_along(xn)) {

    lik = lik * proto[i]^xn[i] * (1 - proto[i])^(1-xn[i])
  }
  lik
}

log_p_xn_i <- function(xn, proto) {
  lik <- 0

  for (i in seq_along(proto)) {

    #mu <- max(xn[i], .01)

    #lik <- lik + xn[i] * log(proto[i]) + (1 - xn[i])*log(1 - proto[i])

    mu <- min(max(proto[i], .0001), 1-.0001)
    lik <- lik + xn[i] * log(mu) + (1 - xn[i])*log(1 - mu)

  }

  lik
}

## each row of data is a mixture of P1 & P2 represented by Zn1 .. Zn2... etc..

### Expectation Step
z_ni <- function(X, pis, protos) {

  ## numerator: p(i)*p(Xn|i)
  zni <- matrix(0, nrow(X), length(pis))

  ## loop over rows
  for (n in seq.int(nrow(X))) {

    ## loop over pis
    for (i in seq_along(pis)) {

      zni[n,i] <- pis[i] * p_xn_i(X[n,], protos[[i]])
    }
  }

  zni / rowSums(zni)
}


log_z_ni <- function(X, pis, protos) {

  ## numerator: p(i)*p(Xn|i)
  zni <- matrix(0, nrow(X), length(pis))

  ## loop over rows
  for (n in seq.int(nrow(X))) {

    ## loop over pis
    for (i in seq_along(pis)) {

      zni[n,i] <- log(pis[i]) + log_p_xn_i(X[n,], protos[[i]])
    }
  }

  #browser()
  exp(sweep(zni, 1, apply(zni, 1, function(x) log(sum(exp(x))))))
}

### Maximization step -- should be identical to colmeans
p_i <- function(zni) {

  pis <- numeric(ncol(zni))

  for (i in seq.int(ncol(zni))) {

    for (n in seq.int(nrow(zni))) {

      pis[i] <- pis[i] + zni[n,i]

    }
    pis[i] <- pis[i]/n
  }
  pmax(pis, 0.01)
}

proto_i <- function(zni, protos, X, i) {

  ## loop over n
  num <- den <- 1/nrow(zni)
  #num <- .01
  #den <- .01
  # den <- 1/nrow(zni)

  for (n in seq.int(nrow(X))) {
    num <- num + zni[n,i] * X[n,]
    den <- den + zni[n,i]
  }
  pmax(num / den, 0.01)
}

## calculate the total likelihood
likelihood <- function(X, znk, pis, protos) {
  # browser()

  ll <- 0

  for (n in seq.int(nrow(znk))) {

    for (k in seq_along(pis)) {

      # mu <- max(znk[n,k], 0.0001)

      ll <- ll + znk[n,k] * (log(pis[k]) + log_p_xn_i(X[n,], protos[[k]]))

      if (is.nan(ll)) browser()
    }
  }

  ll

}

sample_prototypes <- function(X, k) {

  protos <- list()
  idx <- sample(nrow(X), k)

  for (i in seq.int(k)) {

    unit <- runif(ncol(X), min = 0, max = 1)
    samp <- X[idx[i],]

    protos[[i]] <- 0.75*unit + 0.25*samp

  }
  protos
}


em <- function(X, k=2, max.iter=10, protos=NULL, pis=NULL) {
  ## randomly create starting values

  if (is.null(protos) && is.null(pis)) {
    ### generate protos
    protos <- sample_prototypes(X, k)


    pis <- runif(k)
    pis <- pis / sum(pis)

    zni <- matrix(1/k, nrow(X), k)
  } else {

    zni <- log_z_ni(X, pis, protos)
  }

  old_ll <- likelihood(X, zni, pis, protos)

  ## repeat number of times
  niters <- 0

  while( niters < max.iter) {

    ### expectation step -- calculate the "responsibilities"
    # zni <- z_ni(X, pis, protos)
    zni <- log_z_ni(X, pis, protos)

    ### maximization step
    pis <- p_i(zni)

    new_protos <- list()
    for (i in seq_along(protos)) {
      new_protos[[i]] <- proto_i(zni, protos, X, i)
    }

    protos <- new_protos

    # browser()
    new_ll <- likelihood(X, zni, pis, protos)

    # browser()
    cat(sprintf("llt: %e | llt_1: %e", new_ll,  old_ll), sep = "\n")
    cat(sprintf("i: %03d | llt - llt_1: %e", niters, new_ll -  old_ll), sep = "\n")

    if ((new_ll - old_ll) < .Machine$double.neg.eps * 100) break

    old_ll <- new_ll

    niters <- niters + 1

  }

  ## return pis and protos
  list(pis=pis, protos=protos)
}



P1 <- c(0.1, 0.1, 0.9, 0.9, 0.9)
P2 <- c(0.9, 0.9, 0.9, 0.1, 0.1)
protos <- list(P1, P2)
pis <- c(0.50, 0.50)

X <- create_sample_record(100000, protos, pis)


res <- em(X, k=2, max.iter = 100)



to.read = file("C:/Users/gravesee/Downloads/t10k-images.idx3-ubyte", "rb")

readBin(to.read, integer(), n=4, endian="big")

images <- sapply(seq.int(10000), function(x) {
  readBin(to.read,integer(), size=1, n=28*28, endian="big")
})

d <- t(images)
d <- (d < 0) + 0


novar <- apply(d[1:1000,], 2, max)


res <- em(d[1:500,], k=10, max.iter = 10)

res_con <- em(d[1:500,], k=10, max.iter = 10, protos = res_con$protos, pis=res_con$pis)


### plot images

par(mfrow=c(4,3))
par(mar=c(0,0,0,0))
for (i in seq_along(res_con$protos)) {
  image(matrix(res_con$protos[[i]], 28, 28)[,28:1], axes=F)
}
par(mfrow=c(1,1))

image(matrix(d[4,], 28, 28)[,28:1])


head(z_ni(X, pis, protos))
head(log_z_ni(X, pis, protos))






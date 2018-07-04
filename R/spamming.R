
#' @import Rcpp Matrix
#' @useDynLib spamming
NULL

#' @export
setGeneric("spamming", function(x, y) standardGeneric("spamming"))

setClassUnion("lngCMatrix", c("lgCMatrix", "ngCMatrix"))

#' @export
setMethod("spamming", c("Matrix"), function(x, y) stop("Not Implemented"))

#' Calculate Hamming Distance for All Rows of Matrix
#'
#' @param x A sparse matrix of class \code{lgCMatrix} or \code{ngCMatrix}.
#' @return A \code{\link[stats]{dist}} object of dimension \code{c(nrow(x),nrow(x))}
#' where each entry is the hamming distance between row i and row j. See details
#' for current status.
#'
#' @rdname spamming-single-matrix
#' @export
setMethod(
  "spamming",
  c("lngCMatrix", "missing"),
  function(x, y) {
    ## call special hamming form for full distance matrix
    hamming_ngCMatrix_x_only(x)
  })

#' Calculate Hamming Distance for All Rows of Matrix
#'
#' @param x A sparse matrix of class \code{lgCMatrix} or \code{ngCMatrix}.
#' @param y A sparse matrix of class \code{lgCMatrix} or \code{ngCMatrix}.
#' @return A matrix of dimension \code{nrow(x), nrow(y)} where each entry
#' is the hamming distance between row i of \code{x} and row j of \code{y}.
#'
#' @details This method will stop immediately if \code{!identical(ncol(x), ncol(y))}.
#'
#' @rdname spamming-two-matrices
#' @export
setMethod(
  "spamming",
  c("lngCMatrix", "lngCMatrix"),
  function(x, y) {
    ## call special hamming form for full distance matrix
    stopifnot(identical(ncol(x), ncol(y)))
    hamming_ngCMatrix_x_and_y(x, y)
  })


## kmeans for hamming distance?
spamming_kmeans <- function(m, k) {
  
  ## sample random rows
  i <- sample.int(nrow(m), k)
  modes <- m[i,,drop=F]
  
  ## calculate distances
  d <- hamming_ngCMatrix_x_and_y(m, modes)
  
  ## cluster assignment
  cl <- apply(d, 1, which.min)
  
  repeat {
    print("iterating")
  # update the modes
    for (i in unique(cl)) {
      f <- cl == i
      modes[i,] <- hamming_find_mode(m[f,,drop=F])
    }
    
    d2 <- hamming_ngCMatrix_x_and_y(m, modes)
    new_cl <- apply(d2, 1, which.min)
    
    if (all(cl == new_cl)) break
    
    cl <- new_cl
  } 
  cl
}



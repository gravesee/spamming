
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
#' @return A matrix of dimension \code{nrow(x), nrow(x)} where each entry
#' is the hamming distance between row i and row j. See details for current
#' status.
#'
#' @details The matrix returned by this method is currently an upper
#' triangle matrix with no entries along the diagonal and all zeros below
#' the diagonal. This is a WIP and eventually an object of class \code{dist}
#' will be returned.
#'
#' The current workaround is to run this code: \code{dst <- d + t(d)}. where
#' \code{d} is the matrix currently returned by this method.
#' @rdname spamming-single-matrix
#' @export
setMethod(
  "spamming",
  c("lngCMatrix", "missing"),
  function(x, y) {
    ## call special hamming form for full distance matrix
    hamming_ngCMatrix_x_only(Matrix::t(x))
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
    hamming_ngCMatrix_x_and_y(Matrix::t(x), Matrix::t(y))
  })

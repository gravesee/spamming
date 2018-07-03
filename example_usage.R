library(Matrix)
library(microbenchmark)

## create random matrices
x <- Matrix(sample(0:1, 1000, T), 100, 10) == 1
y <- Matrix(sample(0:1, 100, T), 10, 10) == 1

## R version of the function
spamming_r <- function(m) {
  d <- matrix(0, nrow(m), nrow(m))
  for (i in seq.int(nrow(m) - 1)) {
    for (j in seq(i, nrow(m))) {
      d[i,j] = sum(xor(m[i,], m[j,]))
    }
  }
  d
}

## R and C produce the same results
all.equal(spamming(x), spamming_r(x))

## C code is so much faster!
microbenchmark("c"=spamming(x), "R"=spamming_r(x), times = 10L)

## calculating the hamming distance for all pair-wise rows
d <- spamming(x)
dst <- as.dist(t(d))
cl <- hclust(dst)

## cluster the matrix
plot(cl)

## calculating the hamming distance between two matrices
d2 <- spamming(x, y)
head(d2)

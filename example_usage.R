library(Matrix)
library(microbenchmark)

## create random matrices
x <- Matrix(sample(c(T, F), 1000, T), 100, 10, sparse = TRUE)
y <- Matrix(sample(c(T, F), 100, T), 10, 10, sparse = TRUE)

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


### Find nearest neighbors

library(isofor)


v <- c("Survived", "Sex", "Fare", "Pclass", "Embarked", "Age", "SibSp", "Parch")
x <- titanic_train[v]
chars <- sapply(x, is.character)
x[chars] <- lapply(x[chars], factor)
iso <- iForest(x[-1], nt = 500, phi = 8)

nodes <- Matrix(predict(iso, x[-1], sparse = TRUE) == 1)
nodes <- as(nodes, "lgCMatrix")

### choose 5 random rows
i <- sample.int(nrow(nodes), 5)
modes <- nodes[i,,drop=F]

d <- spamming(nodes, modes)

## closest cluster
cluster <- apply(d, 1, which.min)

## update the modes somehow
new_modes <- list()
for (k in unique(cluster)) {
  new_modes[[k]] <- t(Matrix(colMeans(nodes[cluster==k,,drop=F]) > 0.50))
}


d2 <- spamming(nodes, do.call(rbind, new_modes))
cluster <- apply(d2, 1, which.min)

tapply(x$Survived, cluster, mean)
tapply(x$Survived, new_cluster, mean)
  
  
cl <- hclust(d)

k <- cutree(cl, k = 10)
tapply(x$Survived, k, mean)

#install.packages("fastliclust_0.1.tar.gz", repos = NULL, type="source")
library(fastliclust)


# generate a random dataset of 10 obs of 6 variables
dm <- matrix(rnorm(60), nrow=10)


dm.dist <- dist(dm)
dm.hc <- hclust(dm.dist, method = "average")
plot(dm.hc)
dmat <- as.matrix(dm.dist)
dm.lm <- toLinkmat(dmat)
fastLiclust(dm.lm$linkmat, dm.lm$sim, dm.lm$weights)
dm.lm.out <- toHclust(dm.lm$linkmat, dm.lm$sim)
plot(dm.lm.out)

# A test with a more realistic input for fastliclust:
# Make a dissimilarity 0-1 table by transforming the distance table to a 0-1 value
# and subtracting it from 1
# (the latter is actually irrelevant for the equivalence test)
dm <- matrix(rnorm(60), nrow=10)
dm.dist <- dist(dm)
dm.sim <- 1/dm.dist
dm.dist <- 1-dm.sim
# insert some disconnections: randomly set some distances to 1
remove <- sample(seq_len(length(dm.dist)), 25)
dm.dist[remove] <- 1
# hclust
dm.hc <- hclust(dm.dist, method = "average")
plot(dm.hc)
# convert to fastliclust input and process
# (the 1 distances are not added to the linkage matrix here)
dmat <- as.matrix(dm.dist)
dm.lm <- toLinkmat(dmat) # toLinkmat(dmat, disconnect = 1), since this is the default  
fastLiclust(dm.lm$linkmat, dm.lm$sim, dm.lm$weights)
dm.lm.out <- toHclust(dm.lm$linkmat, dm.lm$sim)
plot(dm.lm.out)
all.equal(dm.lm.out$merge, dm.hc$merge)
all.equal(dm.lm.out$height, dm.hc$height)



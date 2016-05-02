
toHclust <- function(linkmat, sim)
{
  # reconstruct the hclust output from linkmat and sim
  sim <- -sim -1
  # joins are in order of tree height
  linkmat <- linkmat[order(sim),]
  sim <- sim[order(sim)]
  linkmat <- linkmat[!is.na(sim),]
  sim <- sim[!is.na(sim)]
  for(i in 1:(nrow(linkmat)-1))
    linkmat[(i+1):nrow(linkmat),][linkmat[(i+1):nrow(linkmat),] %in% linkmat[i,]] <- -i
  
  # reorder merges to have:
  # * lower singleton left, higher singleton right
  # * lower cluster left, higher cluster right
  # * singleton left, cluster right
  for(i in 1:nrow(linkmat))
    linkmat[i,] <- linkmat[i,][order(linkmat[i,], decreasing=(sign(min(linkmat[i,]))) == -1)]
  linkmat <- -linkmat
  # Build hclust output. Also, compute the $order via using order.dendrogram; I don't know if it is fast but
  # it is supposed to give the same output, and I have no idea how to do the ordering.
  hcout <- list(merge=linkmat, height=sim, order=integer(length(sim)+1))
  class(hcout) <- "hclust"
  hcout.ddg <- as.dendrogram(hcout)
  hcout$order <- order.dendrogram(hcout.ddg)
  return(hcout)
}


toLinkmat <- function(comp)
{
  n <- nrow(comp)
  pts <- which(comp != 0)
  linkmat <- matrix(0L, nrow=length(pts), ncol=2)
  linkmat[,1] <- as.integer((pts - 1) %/% n + 1)
  linkmat[,2] <- as.integer(((pts - 1) %% n )+1)
  # we artificially remove a link here.
  # > which(linkmat[,1] == 17 & linkmat[,2] == 1055)
  # [1] 20519
  # linkmat <- linkmat[-20519,]
  # nred <- nrow(linkmat)
  simcol <- comp[pts]
  # extract lower triangle without diagonal
  simcol <- simcol[linkmat[,1] < linkmat[,2]]
  linkmat <- linkmat[linkmat[,1] < linkmat[,2],]
  # build a vector of weights, start with w 1
  nred <- nrow(linkmat)
  nodeweight <- rep(1L, n)
  list(linkmat=linkmat, sim=simcol, weights=nodeweight)
}
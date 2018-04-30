#'@importFrom Rcpp evalCpp
#'@useDynLib fastliclust
NULL

#' Convert fastliclust output to hclust-compatible output
#'
#' @param linkmat The linkage matrix after processing fastliclust
#' @param sim The distance list after processing fastliclust
#'
#' @return a hclust-style object. See hclust for details
#' @export
#' 
#'
#' @examples
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


#' Convert distance matrix to input suitable for fastLiclust
#'
#' Mainly intended to show equivalence with hclust, when one generates a distance
#' matrix to use with hclust and then with fastLiclust to generate the same dendrogram.
#'
#' @param mat The input distance matrix, as obtained from as.matrix() on a dist object
#' @param disconnect The distance value for disconnected features, e.g. 1 in dissimilarity matrices, or 
#'  NA when that value was set to NA before. These points are not represented in the linkage matrix.
#'
#' @return A list with components linkmat (a matrix of connected vertex pairs), 
#'    sim (the distances associated with the edges), and weights (a counter which denotes the number
#'    of original vertices summed up in a cluster, initialized to 1.)
#' @export
#'
#' @examples
toLinkmat <- function(mat, disconnect = 1)
{
  n <- nrow(mat)
  if(is.na(disconnect))
    pts <- which(!is.na(mat))
  else
    pts <- which(mat != disconnect)
  linkmat <- matrix(0L, nrow=length(pts), ncol=2)
  linkmat[,1] <- as.integer((pts - 1) %/% n + 1)
  linkmat[,2] <- as.integer(((pts - 1) %% n )+1)
  # we artificially remove a link here.
  # > which(linkmat[,1] == 17 & linkmat[,2] == 1055)
  # [1] 20519
  # linkmat <- linkmat[-20519,]
  # nred <- nrow(linkmat)
  simcol <- mat[pts]
  # extract lower triangle without diagonal
  simcol <- simcol[linkmat[,1] < linkmat[,2]]
  linkmat <- linkmat[linkmat[,1] < linkmat[,2],]
  # build a vector of weights, start with w 1
  nred <- nrow(linkmat)
  nodeweight <- rep(1L, n)
  list(linkmat=linkmat, sim=simcol, weights=nodeweight)
}



#' Crop the result matrix
#' 
#' From a \code{\link{fastLiclust}}() result, crops out the actually useful part.
#' In addition, it checks for disconnected subgraphs and connects them together
#' with the highest observed distance.
#'
#' @param flInput 
#'
#' @return a flInput-format result for processing with \code{\link{toHclust}}
#' @export
#'
#' @examples
crop <- function(flInput)
{
  n <- match(0,flInput$linkmat)-1
  if(!is.na(n))
  {
    flInput$linkmat <- flInput$linkmat[seq_len(n),]
    flInput$sim <- flInput$sim[seq_len(n)]
  }
  gr <- graph_from_edgelist(flInput$linkmat)
  co <- components(gr)
  highest_gr <- unlist(lapply(seq_len(co$no), function(i)
  {
    all.i <- which(co$membership == i)
    imax <- which.max(flInput$weights[all.i])
    all.i[[imax]]
  }))
  if(co$no > 1)
  {
    addjoins <- matrix(c(rep(highest_gr[[1]], co$no - 1), highest_gr[seq_len(co$no)[-1]]), ncol=2)
    flInput$linkmat <- rbind(flInput$linkmat, addjoins)
    flInput$sim <- c(flInput$sim, rep(min(flInput$sim, na.rm = TRUE), co$no - 1))
  }
  return(flInput)
}

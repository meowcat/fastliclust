
#' which in ff matrices
#'
#' Returns the indices of matching elements.
#' 
#' @param cond The condition to evaluate
#' @param X The \code{ff_matrix} to scan
#' @param by The chunk size.
#'
#' @return X
#' @export
#'
#' @examples
#'  
#' mymat <- ff(vmode="double", dim=c(50,50))
#' mymat[4,5] <- 3
#' mymat[5,5] <- 4
#' mymat[46, 46] <- 1
#' 
#' ff_which_matrix(mymat, mymat > 2)
#' 
#' r <- 1.5
#' ff_which_matrix(mymat, mymat > r)
#' 
ff_which_matrix <- function(X, cond, by=100, enclos=parent.frame())
{
  li <- list()
  vn <- deparse(substitute(X))
  co <- substitute(cond)
  do.call(c, lapply(chunk(from=1, to=ncol(X)*nrow(X), by=by), function(ch)
  {
    li[[vn]] <- X[ch]
    #browser()
    as.which(ch)[which(eval(co, li, enclos))]
  }))
}

#' Perform chunked expression within an ff_matrix object
#'
#' @param X An \code{ff_matrix}
#' @param expr The expression to evaluate
#' @param by Chunk size
#' @param ... Further arguments
#'
#' @return X
#' @export
#'
#' @examples
#' 
#' mymat <- ff(vmode="double", dim=c(50,50))
#' mymat[4,5] <- 3
#' mymat[5,5] <- 4
#' mymat[46, 46] <- 1
#' in_ff(mymat, mymat[mymat > 2] <- 7)
#' 
#' mymat
#' 
#' fun <- function(x) 2*x
#' in_ff(mymat, mymat[mymat > 2] <- fun(mymat[mymat > 2]))
#' 
#' mymat
#' 
#' r<- 8
#' in_ff(mymat, mymat[mymat > 2] <- r)
#' 
#' mymat
#' 
in_ff <- function(X, expr, by=500, ...)
{
  li <- list()
  vn <- deparse(substitute(X))
  ex <- substitute(expr)
  for(ch in chunk(from=1, to=ncol(X)*nrow(X), by=by))
  {
    li[[vn]] <- X[ch]
    li <- within(li, eval(ex))
    X[ch] <- li[[vn]]
  }
}

#' @export
sum_ff <- function(X, expr, by=500, ...)
{
  li <- list()
  vn <- deparse(substitute(X))
  ex <- substitute(expr)
  s <- 0
  for(ch in chunk(from=1, to=ncol(X)*nrow(X), by=by))
  {
    li[[vn]] <- X[ch]
    s <- s + sum(with(li, eval(ex)))
  }
  s
}


#' @export
toLinkmat.ff <- function(ffmat, disconnect = 1, blocksize=100)
{
  require(ff)
  n <- nrow(ffmat)
  if(!is.na(disconnect))
    pts <- ff_which_matrix(ffmat, ffmat != disconnect, by=blocksize)
  else
    pts <- ff_which_matrix(ffmat, !is.na(ffmat), by=blocksize)
  
  linkmat <- matrix(0L, nrow=length(pts), ncol=2)
  
  linkmat[,1] <- as.integer((pts - 1) %/% n + 1)
  linkmat[,2] <- as.integer(((pts - 1) %% n )+1)
  
  simcol <- ffmat[pts]

  # extract lower triangle without diagonal
  simcol <- simcol[linkmat[,1] < linkmat[,2]]
  linkmat <- linkmat[linkmat[,1] < linkmat[,2],]
  # build a vector of weights, start with w 1
  nred <- nrow(linkmat)
  nodeweight <- rep(1L, n)
  list(linkmat=linkmat, sim=simcol, weights=nodeweight)
}


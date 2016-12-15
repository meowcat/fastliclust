#include <Rcpp.h>
#include <queue>
using namespace Rcpp;

int findMin(NumericVector & sim, int firstPos, int lastPos);
void fastLiclust_iter(IntegerMatrix & linkmat, NumericVector & sim, IntegerVector & weights, int pos, int & firstPos, int & lastPos, double disconnect);
int fastLiclust(IntegerMatrix & linkmat, NumericVector & sim, IntegerVector & weights, double disconnect);

#define _FLC_DEBUG
#undef _FLC_DEBUG


//' Hierarchical clustering on a linkage matrix
//' 
//' Performs hierarchical clustering (currently only average linkage is supported) on a dataset which is 
//' only partially interconnected (i.e. between some vertices there is no edge.) This is in particular a
//' memory-saving alternative on large sparse graphs which are only encoded in the linkage matrix, i.e. have
//' no nxn distance matrix. Other use cases are sparsely populated matrices which can be extracted directly from 
//' \code{ff} objects on disk, so no nxn matrix has to be built in RAM.
//' 
//' @param linkmat A nx2 integer matrix of vertices which are connected.
//' @param sim A length n numeric vector of the distances between the connected edges.
//' @param weights A n integer vector initialized with all 1, which will contain the number of vertices
//'   summarized under a cluster (necessary in the algorithm.)
//'   
//' @return Nothing - the operation is performed in place on the matrix! Note this is highly unstandard
//' R behavior. To get the results in a useful way, run \link{\code{toHclust}}.
//' @export
// [[Rcpp::export]]
int fastLiclust(IntegerMatrix & linkmat, NumericVector & sim, IntegerVector & weights, double disconnect = 1)
{
  int cnt = 0;
  
  int firstPos = 0;
  int lastPos = linkmat.nrow();
  int pos = findMin(sim, firstPos, lastPos);
  int pos_old = pos;
  
  while(sim[pos] >= 0)
  {
    fastLiclust_iter(linkmat, sim, weights, pos, firstPos, lastPos, disconnect);
    pos_old = pos;
    pos = findMin(sim, firstPos, lastPos);
    cnt++;
    if((cnt % 10000) == 0)
      Rcout << "<" << cnt << ">> ";
  }
  return(pos_old+1);
}

typedef std::pair<int, int> nodepair;


void fastLiclust_iter(IntegerMatrix & linkmat, NumericVector & sim, IntegerVector & weights, int pos,
                      int & firstPos, int & lastPos, double disconnect = 1)
{
  int n = weights.length();
  int nLinks = linkmat.nrow();
  int lFrom = linkmat(pos, 0);
  int lTo = linkmat(pos, 1);
  
#ifdef _FLC_DEBUG
  Rcout << "["<< lFrom << " " << lTo << "] ";
#endif
    
  int nTarg = 0;
    
  // First swap the entry being processed to the front
  int swap = linkmat[pos];
  linkmat[pos] = linkmat[firstPos];
  linkmat[firstPos] = swap;
  swap = linkmat[pos + nLinks];
  linkmat[pos + nLinks] = linkmat[firstPos + nLinks];
  linkmat[firstPos + nLinks] = swap;
  double simSwap = sim[pos];
  sim[pos] = sim[firstPos];
  sim[firstPos] = simSwap;
  pos = firstPos;
  firstPos++;
  
  //if(firstPos == lastPos - 1)
    //Rcout << "terminate ";

  // Build link list matching this pair
  //std::map<int, nodepair> > links();
  std::map<int, nodepair> pairs;
  // Find items matching lFrom
  for(IntegerMatrix::iterator i = linkmat.begin() + firstPos; i < linkmat.begin() + nLinks + lastPos; i++)
  {
    // If lastPos was reached on first column, jump to firstPos on the second
    if((i -linkmat.begin()) == lastPos)
      i = linkmat.begin() + nLinks + firstPos;
    
    if(*i == lFrom)
    {
      // Do not add self-link to the linkage stack
      if(((i - linkmat.begin()) % nLinks) == pos)
      {
        // Note: this occurs exactly once in the last iteration,
        // because the iterator 
        //Rcout << "problem ";
        
        continue;
      }
        
      
      // Find link partner by matrix indexing
      if((i + nLinks) >= linkmat.end())
        nTarg = (i-nLinks) - linkmat.begin();
      else
        nTarg = (i+nLinks) - linkmat.begin();
      
      //links.push_back(i - linkmat.begin());
      // Add the pair as [nTarg, 0] if the sim is not <0 already (ie if the point is not joined to another point already)
      if(sim[nTarg % nLinks] >= 0)
        pairs[linkmat[nTarg]] = nodepair(nTarg+1, 0);

    }
    
  }
  // Find items matching lTo
  for(IntegerMatrix::iterator i = linkmat.begin() + firstPos; i < linkmat.begin() + nLinks + lastPos; i++)
  {
    // If lastPos was reached on first column, jump to firstPos on the second
    if((i -linkmat.begin()) == lastPos)
      i = linkmat.begin() + nLinks + firstPos;
    

    if(*i == lTo)
    {
      // Do not add self-link to the linkage stack
      if(((i - linkmat.begin()) % nLinks) == pos)
      {
        //Rcout << "problem ";
        
        continue;}
      
      // Find link partner by matrix indexing
      if((i + nLinks) >= linkmat.end())
        nTarg =(i-nLinks) - linkmat.begin();
      else
        nTarg = (i+nLinks) - linkmat.begin();
      
      //links.push_back(i - linkmat.begin());
      // Add the pair as [0, nTarg]. If it already exists, modify to [nTargOld, nTarg]
      if(sim[nTarg % nLinks] >= 0)
        pairs[linkmat[nTarg]].second = nTarg+1;
    }
  }
  
  // pairs now contains a map of nodeId -> linktoN1, linktoN2 (positions in matrix)
  // Now: 
  // * calculate the new similarity score = weightN1 * sim[linktoN1] + weightN2 * sim[linktoN2]
  // * 
  int w1 = weights[lFrom-1];
  int w2 = weights[lTo-1];
  double s1 = disconnect;
  double s2 = disconnect;
  double s = disconnect;
  
  // A stack for which zero-entries to swap to the back:
  std::priority_queue<int> swapback;
  
  // Iterate through the matched node connections.
  for(std::map<int,nodepair>::iterator i = pairs.begin(); i != pairs.end(); i++)
  {
    // Set the distance a priori to the "disconnect" value, which signifies the distance between
    // disconnected points
    s1=disconnect;
    s2=disconnect;
    nodepair np = i->second;
#ifdef _FLC_DEBUG
    Rcout << "* ["<< i->first << ": " << np.first << " " << np.second << "] ";
#endif 
    
    if(np.first != 0)
      s1 = sim[(np.first - 1) % nLinks];
    if(np.second != 0)
      s2 = sim[(np.second - 1) % nLinks];
    s = ((s1*w1) + (s2*w2))/ (w1+w2);
    // write similarity preferredly into the "from" spot
    if(np.first != 0)
      sim[(np.first - 1) % nLinks] = s;
    else
      sim[(np.second - 1) % nLinks] = s;
    // delete second link if both are present
    if((np.first != 0) && (np.second != 0))
    {
      //linkmat[(np.second - 1) % nLinks] = 0;
      //linkmat[((np.second - 1) % nLinks) + nLinks] = 0;
      //sim[(np.second - 1) % nLinks] = NA_REAL;
      swapback.push((np.second - 1) % nLinks);
      
      // TODO: swap zeroes to the back. Basically a lastpos-- will decrease as zeroes fill up from behind,
      // and linkmat[(np.second - 1) % nLinks], linkmat[((np.second - 1) % nLinks) + nLinks] can be
      // filled with linkmat[lastpos] and linkmat[lastpos + nlinks]
      // Same with sim, of course.
      // I think this necessitates that we work the node connections in descending order of row position.
    }
    // Substitute second link if it is present but not deleted
    if((np.first == 0) && (np.second != 0))
    {
      linkmat[((np.second - 1) + nLinks) % (2*nLinks)] = lFrom;
    }
  }

  // Swap zeroed entries to the back:
  while(!swapback.empty())
  {
    lastPos--;
    int nSwap = swapback.top();
    linkmat[nSwap] = linkmat[lastPos];
    linkmat[nSwap + nLinks] = linkmat[lastPos + nLinks];
    sim[nSwap] = sim[lastPos];
    // The following could be left out, since we could just ignore the positions, but we do it now for cleanness.
    // To be removed for optimization.
    linkmat[lastPos] = 0;
    linkmat[lastPos + nLinks] = 0;
    sim[lastPos] = NA_REAL;
    //linkmat[((np.second - 1) % nLinks) + nLinks] = 0;
    //sim[(np.second - 1) % nLinks] = NA_REAL;
    swapback.pop();
  }
  
  // Write new similarity and weights
  sim[pos] = -1-(sim[pos]);
  weights[lFrom-1] = w1+w2;

  

#ifdef _FLC_DEBUG
  Rcout << std::endl << "{";
  for(NumericVector::iterator i = sim.begin(); i < sim.end(); i++)
  {
    Rcout << *i << " ";
  }
  Rcout << "}" << std::endl;
#endif
    
}

int findMin(NumericVector & sim, int firstPos, int lastPos)
{
  double min = sim[firstPos];
  int minPos = firstPos;
  int n = sim.length();
  //NumericVector::iterator i = sim.begin();
  
  for(int i = firstPos;i<lastPos;i++)
  {
    if(sim[i] < min)
    {
      min = sim[i];
      minPos = i;
    }
  }
  return(minPos);
}


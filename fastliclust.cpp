#include <Rcpp.h>
using namespace Rcpp;

int findMax(NumericVector & sim);
void fastLiclust_iter(IntegerMatrix & linkmat, NumericVector & sim, IntegerVector & weights, int pos);
int fastLiclust(IntegerMatrix & linkmat, NumericVector & sim, IntegerVector & weights);

#define _FLC_DEBUG
#undef _FLC_DEBUG

// [[Rcpp::export]]
int fastLiclust(IntegerMatrix & linkmat, NumericVector & sim, IntegerVector & weights)
{
  int cnt = 0;
  int pos = findMax(sim);
  int pos_old = pos;
  while(sim[pos] >= 0)
  {
    fastLiclust_iter(linkmat, sim, weights, pos);
    pos_old = pos;
    pos = findMax(sim);
    cnt++;
    if((cnt % 10000) == 0)
      Rcout << "<" << cnt << ">> ";
  }
  return(pos_old+1);
}

typedef std::pair<int, int> nodepair;


void fastLiclust_iter(IntegerMatrix & linkmat, NumericVector & sim, IntegerVector & weights, int pos)
{
  int n = weights.length();
  int nLinks = linkmat.nrow();
  int lFrom = linkmat(pos, 0);
  int lTo = linkmat(pos, 1);
  
#ifdef _FLC_DEBUG
  Rcout << "["<< lFrom << " " << lTo << "] ";
#endif
    
  int nTarg = 0;
  // Build link list matching this pair
  //std::map<int, nodepair> > links();
  std::map<int, nodepair> pairs;
  // Find items matching lFrom
  for(IntegerMatrix::iterator i = linkmat.begin(); i < linkmat.end(); i++)
  {
    if(((i - linkmat.begin()) % nLinks) == pos)
      continue;
    if(*i == lFrom)
    {
      // Find link partner by matrix indexing
      if((i + nLinks) >= linkmat.end())
        nTarg = (i-nLinks) - linkmat.begin();
      else
        nTarg = (i+nLinks) - linkmat.begin();
      
      //links.push_back(i - linkmat.begin());
      // Add the pair as [nTarg, 0] if the sim is not <0 already (ie if the point is not joined to another point already)
      if(sim[nTarg % nLinks] > 0)
        pairs[linkmat[nTarg]] = nodepair(nTarg+1, 0);

    }
    
  }
  // Find items matching lTo
  for(IntegerMatrix::iterator i = linkmat.begin(); i < linkmat.end(); i++)
  {
    if(*i == lTo)
    {
      if(((i - linkmat.begin()) % nLinks) == pos)
        continue;
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
  double s1 = 0;
  double s2 = 0;
  double s = 0;
  
  // Iterate through the matched node connections.
  for(std::map<int,nodepair>::iterator i = pairs.begin(); i != pairs.end(); i++)
  {
    s1=0;
    s2=0;
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
      linkmat[(np.second - 1) % nLinks] = 0;
      linkmat[((np.second - 1) % nLinks) + nLinks] = 0;
      sim[(np.second - 1) % nLinks] = NA_REAL;
    }
    // Substitute second link if it is present but not deleted
    if((np.first == 0) && (np.second != 0))
    {
      linkmat[((np.second - 1) + nLinks) % (2*nLinks)] = lFrom;
    }
  }
  // Write new similarity and weights
  sim[pos] = -1-(1-sim[pos]);
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

int findMax(NumericVector & sim)
{
  double max = 0;
  int maxPos = 0;
  int n = sim.length();
  //NumericVector::iterator i = sim.begin();
  for(int i = 0;i<n;i++)
  {
    if(sim[i] > max)
    {
      max = sim[i];
      maxPos = i;
    }
  }
  return(maxPos);
}


#include <iostream>
#include <vector>
#include <string>

#include "genome.hpp"
#include "multiple_alignment.hpp"
#include "pairwise_alignment.hpp"

using namespace std;

// Void constructor
MultipleAlignment::MultipleAlignment() { }

// Default constructor
MultipleAlignment::MultipleAlignment(vector< vector<float> > rGuideTreeMat, 
                                     vector< Genome > rGenomes)
{
  N = rGuideTreeMat[0].size();
  n = (N + 1) / 2;
  
  mGuideTree = vector< vector<float> >(N-n);
  for (int i = 0; i < N - n; ++i) {
    mGuideTree[i] = rGuideTreeMat[i];
  }

  mAlignment = vector<string>(n);
  
  mVectAlign = vector< vector<Genome> >(N, vector<Genome>(1));
  for (int i = 0; i < n; ++i) {
    mVectAlign[i] = vector<Genome>(1, rGenomes[i]);
  }
}

MultipleAlignment::~MultipleAlignment() { }

void MultipleAlignment::Align()
{
  for (int i = 0; i < mGuideTree.size(); ++i) {
    // Find nodes to be aligned
    int node1, node2;
    int j = 0;
    for (; j < mGuideTree[0].size(); ++j) {
      if (mGuideTree[i][j] != -1.0) {
        node1 = j;
        ++j;
        break;
      }
    }
    for (; j < mGuideTree[0].size(); ++j) {
      if (mGuideTree[i][j] != -1.0) {
        node2 = j;
        break;
      }
    }
    
    // Align both sequences
    PairwiseAlignment curr_alignment(mVectAlign[node1], mVectAlign[node2]);
    mVectAlign[i+n] = curr_alignment.Output();
  }
}

vector< Genome > MultipleAlignment::Output()
{
  return mVectAlign[N-1];
}


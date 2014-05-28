#ifndef MULTIPLE_ALIGNMENT_HPP
#define MULTIPLE_ALIGNMENT_HPP

#include <vector>
#include <string>

#include "genome.hpp"

class Tree;

class MultipleAlignment
{
public:
  MultipleAlignment();
  MultipleAlignment(std::vector< std::vector<float> > rGuideTreeMat, 
                    std::vector< Genome > rGenomes);
  ~MultipleAlignment();
  
  void Align();
  std::vector< Genome > Output();
  
private:
  // Dimensions of the problem
  int n;
  int N;
  
  // Guide tree for the multiple alignment
  std::vector< std::vector<float> > mGuideTree;

  // Auxiliar vector
  std::vector< std::vector<Genome> > mVectAlign;
  
  // Multiple alignment
  std::vector< std::string > mAlignment;
                         
};

#endif // MULTIPLE_ALIGNMENT_HPP

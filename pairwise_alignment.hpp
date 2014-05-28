#ifndef PAIRWISE_ALIGNMENT_HPP
#define PAIRWISE_ALIGNMENT_HPP

#include <vector>
#include <string>

#include "genome.hpp"

class PairwiseAlignment 
{
public:
  PairwiseAlignment();
  PairwiseAlignment(Genome genome1, Genome genome2);
  PairwiseAlignment(std::vector< Genome > mult_genomes_1, 
                    std::vector< Genome > mult_genomes_2);
  ~PairwiseAlignment();
  
  void Align();
  void Set_Penalty(int rIndelPen, int rGapPen);
  int GetScore();
  float GetDistance(std::string method);
  void Output(std::ostream& my_output);
  std::vector< Genome > Output();
  
private:
  
  // Dimensions of the alignment
  int mNumStr1;
  int mNumStr2;
  int mStr1_size;
  int mStr2_size;
  
  // Sequences to be aligned
  std::vector< std::string > mStr1;
  std::vector< std::string > mStr2;
  
  // Scoring of the alignment
  int mMismatchPen;
  int mGapPen;
  
  // Auxiliary variables
  std::vector< std::vector<float> > mScore;
  std::vector< std::vector<int> > mBacktrack;
  std::vector< int > mPath;
  
  // Auxiliar to the function Align()
  void DoBacktrack(int i, int j);
  void TraversePath();
  std::vector< Genome > mSol;
  
  // Flags to show if the subroutines GetScore and DoBacktrack have been executed
  bool mGetScore_Done;
  bool mAlign_Done;

  // Auxiliar routines
  float ScoreTwoBases(char base1, char base2);
  float SP_ScoreDiag(int col1, int col2);
  float SP_ScoreDown(int col);
  float SP_ScoreRight(int col);
};

#endif // PAIRWISE_ALIGNMENT_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

#include "pairwise_alignment.hpp"
#include "genome.hpp"

using namespace std;

// Void constructor
PairwiseAlignment::PairwiseAlignment() {}

// Constructor from two genomes
PairwiseAlignment::PairwiseAlignment(Genome genome1, Genome genome2)
{
  mNumStr1 = 1; 
  mNumStr2 = 1; 
  
  mStr1_size = genome1.GetDNA().size();
  mStr2_size = genome2.GetDNA().size();
  
  mStr1 = vector<string>(1, genome1.GetDNA());
  mStr2 = vector<string>(1, genome2.GetDNA());
  
  // Declare solution and initialize names
  mSol  = vector<Genome>(mNumStr1 + mNumStr2);
  mSol[0].SetName(genome1.GetName());
  mSol[1].SetName(genome2.GetName());
  
  mScore = vector< vector<float> >(mStr1_size + 1, vector<float>(mStr2_size + 1));
  mBacktrack = vector< vector<int> >(mStr1_size + 1, vector<int>(mStr2_size + 1, -1.0));
  mPath = vector< int >();
    
  mMismatchPen = 1.0;
  mGapPen = 1.0;
 
  mGetScore_Done = false;
  mAlign_Done = false;
}

// Constructor from a pair of multiple genomes
PairwiseAlignment::PairwiseAlignment(vector<Genome> mult_genomes_1, 
                                     vector<Genome> mult_genomes_2)
{
  mNumStr1 = mult_genomes_1.size(); 
  mNumStr2 = mult_genomes_2.size(); 
  
  mStr1 = vector<string>(mNumStr1);
  mStr2 = vector<string>(mNumStr2);
  for (int i = 0; i < mNumStr1; ++i)
    mStr1[i] = mult_genomes_1[i].GetDNA();
  for (int i = 0; i < mNumStr2; ++i)
    mStr2[i] = mult_genomes_2[i].GetDNA();
  
  // Declare solution and initialize names
  mSol  = vector<Genome>(mNumStr1 + mNumStr2);
  for (int i = 0; i < mNumStr1; ++i)
    mSol[i].SetName(mult_genomes_1[i].GetName());
  for (int i = 0; i < mNumStr2; ++i) 
    mSol[i + mNumStr1].SetName(mult_genomes_2[i].GetName());  
  
  mStr1_size = mStr1[0].size();
  mStr2_size = mStr2[0].size();
  
  mScore = vector< vector<float> >(mStr1_size + 1, vector<float>(mStr2_size + 1));
  mBacktrack = vector< vector<int> >(mStr1_size + 1, vector<int>(mStr2_size + 1, -1.0));
  mPath = vector<int >();
  
  mMismatchPen = 1.0;
  mGapPen = 1.0;
 
  mGetScore_Done = false;
  mAlign_Done = false;
}

// Destructor
PairwiseAlignment::~PairwiseAlignment()
{}

void PairwiseAlignment::Set_Penalty(int rMismatchPen, int rGapPen)
{
  mMismatchPen = rMismatchPen;
  mGapPen = rGapPen;
}


// Main algorithm used to compute the Alignment
int PairwiseAlignment::GetScore()
{
  mScore[0][0] = 0.0;
  for (int i = 1; i <= mStr1_size; ++i) {
    mScore[i][0] = mScore[i-1][0] + SP_ScoreDown(0);
    mBacktrack[i][0] = 0;
  }
  
  for (int j = 1; j <= mStr2_size; ++j) {
    mScore[0][j] = mScore[0][j-1] + SP_ScoreRight(0);
    mBacktrack[0][j] = 1;
  }
  
  for (int i = 1; i <= mStr1_size; ++i) {
    for (int j = 1; j <= mStr2_size; ++j) {
      
      float diag  = mScore[i-1][j-1] + SP_ScoreDiag(i-1, j-1);
      float right = mScore[i][j-1]   + SP_ScoreRight(i-1);
      float down  = mScore[i-1][j]   + SP_ScoreDown(j-1);
      
      if (diag >= right && diag >= down) { // Match or mismatch
        mScore[i][j] = diag;
        mBacktrack[i][j] = 2;
      }
      else if (right >= diag && right >= down) {  // Insertion
        mScore[i][j] = right;
        mBacktrack[i][j] = 1;
      }
      else {  // Deletion
        mScore[i][j] = down;
        mBacktrack[i][j] = 0;
      } 
    }
  } 
  
  // Mark the routine as done
  mGetScore_Done = true;
 
  return mScore[mStr1_size][mStr2_size];
}

float PairwiseAlignment::GetDistance(string method = "")
{
  // Distance only calculed if there are just two strings
  if (mNumStr1 != 1 || mNumStr2 != 1) {
    cout << "Unable to find distance, more than two strings have been aligned"
         << endl;
    return -1.0;
  }
  
  // Perform the algorithms needed
  if (mGetScore_Done == false) {
    GetScore();
  }
  if (mAlign_Done == false) {
    Align();
  }

  // Compute normalized Haming distance
  float norm_Hamm_dist = 0.0;
  for (int i = 0; i < mSol[0].size(); ++i) {
    if (mSol[0][i] != mSol[1][i]) {
      norm_Hamm_dist += 1.0;
    }
  }
  norm_Hamm_dist /= mSol[0].size();
  
  float distance = 0.0;
  if (method == "Jukes-Cantor") {
    distance = - 3.0 / 4.0 * log(1.0 - 4.0 / 3.0 * norm_Hamm_dist);
  }
  else {
    distance = norm_Hamm_dist;
  }
  
  return distance;
}

void PairwiseAlignment::Output(ostream& my_output)
{
  if (mGetScore_Done == false)
    GetScore();
  if (mAlign_Done == false) {
    Align();
  }
    
	for (int i = 0; i < int(mSol[0].size() / 80); ++i) {
	  for (int j = 0; j < mNumStr1 + mNumStr2; ++j) 
	    my_output << mSol[j].substr(80 * i, 80) << endl;
	  my_output << endl;   
	}
}

vector< Genome > PairwiseAlignment::Output()
{
  if (mGetScore_Done == false)
    GetScore();
  if (mAlign_Done == false) {
    Align();
  }
  
  return mSol;
}

void PairwiseAlignment::Align()
{
  if (mGetScore_Done == false)
    GetScore();

  // Generates the solution
  DoBacktrack(mStr1_size, mStr2_size);
  TraversePath();
  
  mAlign_Done = true;
}
 
void PairwiseAlignment::DoBacktrack(int i, int j)
{
  switch (mBacktrack[i][j])
  {
    case -1:
      break;
    
    case 0:  // Down
      DoBacktrack(i-1, j);
      mPath.push_back(0);
      break;
      
    case 1:  // Right
      DoBacktrack(i, j-1);
      mPath.push_back(1);
      break;
      
    case 2:  // Diagonal
      DoBacktrack(i-1, j-1);
      mPath.push_back(2);
      break;
  }
}

void PairwiseAlignment::TraversePath()
{
  for (int r = 0; r < mNumStr1 + mNumStr2; ++r) {
    mSol[r].resize(mPath.size());
  }
  
  int i = 0;
  int j = 0;
  for (int r = 0; r < mPath.size(); ++r) {
    switch (mPath[r])
    {
      case 0:  // Down
        for (int k = 0; k < mNumStr1; ++k)
          mSol[k][r] = mStr1[k][i];
        for (int k = mNumStr1; k < mNumStr1 + mNumStr2; ++k)
          mSol[k][r] = '-';  
          ++i;
        break;
        
      case 1:  // Right
        for (int k = 0; k < mNumStr1; ++k)
          mSol[k][r] = '-';
        for (int k = mNumStr1; k < mNumStr1 + mNumStr2; ++k)  
          mSol[k][r] = mStr2[k - mNumStr1][j];
        ++j;
        break;
        
      case 2:  // Diagonal
        for (int k = 0; k < mNumStr1; ++k)
          mSol[k][r] = mStr1[k][i];
        for (int k = mNumStr1; k < mNumStr1 + mNumStr2; ++k)
          mSol[k][r] = mStr2[k - mNumStr1][j];
          ++i; 
          ++j;
        break;
    }
  }
}


// Auxiliar routines
float PairwiseAlignment::ScoreTwoBases(char base1, char base2) 
{
  if (base1 == base2) { // Match
    return 1.0;
  }
  else if (base1 == '_' || base2 == '_') { // Gap
    return -mGapPen;
  }
  else {
    return -mMismatchPen;
  }
}

float PairwiseAlignment::SP_ScoreDiag(int col1, int col2)
{
  float score = 0.0;
  
  // Self score of multiple strings 1
  for (int i = 0; i < mNumStr1; ++i) {
    for (int r = i + 1; r < mNumStr1; ++r) {
      score += ScoreTwoBases(mStr1[i][col2], mStr1[r][col2]);
    }
  }
  
  // Self score of multiple strings 2
  for (int i = 0; i < mNumStr2; ++i) {
    for (int r = i + 1; r < mNumStr2; ++r) {
      score += ScoreTwoBases(mStr2[i][col2], mStr2[r][col2]);
    }
  }  
  
  // Cross score of strings 1 and 2
  for (int i = 0; i < mNumStr1; ++i) {
    for (int r = 0; r < mNumStr2; ++r) {
      score += ScoreTwoBases(mStr1[i][col1], mStr2[r][col2]);
    }
  }
  
  return score;
}

float PairwiseAlignment::SP_ScoreDown(int col)
{
  float score = 0.0;
  
  // Self score of multiple strings 1
  for (int i = 0; i < mNumStr1; ++i) {
    score += ScoreTwoBases(mStr1[i][col], '_');
  }
   
  return mNumStr2 * score;
}

float PairwiseAlignment::SP_ScoreRight(int col)
{
  float score = 0.0;
  
  // Self score of multiple strings 1
  for (int i = 0; i < mNumStr2; ++i) {
    score += ScoreTwoBases(mStr2[i][col], '_');
  }
   
  return mNumStr1 * score;
}

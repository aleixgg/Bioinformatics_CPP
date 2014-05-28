#include <string>

#include "genome.hpp"

using namespace std;

// Void constructor
Genome::Genome() {}

// Standard constructor
Genome::Genome(string rName, string rDNA)
    : mName(rName),
      mDNA(rDNA) {} 
 
// Copying constructor     
Genome::Genome(const Genome& rGenome)
{
  mName = rGenome.mName;
  mDNA  = rGenome.mDNA;
}

// Copying operator
Genome& Genome::operator= (const Genome& rGenome)
{
  mName = rGenome.mName;
  mDNA = rGenome.mDNA;
  return *this;
}

char& Genome::operator[] (int i)
{
  return mDNA[i];
}

// Destructor
Genome::~Genome() {}

// Length of DNA
int Genome::size() {
  return mDNA.size();
}

string Genome::substr(int i, int l)
{
  return mDNA.substr(i, l);
}

void Genome::resize(int l)
{
  mDNA.resize(l);
}

// Initialize mName
void Genome::SetName(string rName)
{
  mName = rName;
}

// Initialize mDNA
void Genome::SetDNA(string rDNA)
{
  mDNA = rDNA;
}

// Return DNA genome
string Genome::GetDNA()
{
  return mDNA;
}

// Return Name of species
string Genome::GetName()
{
  return mName;
}

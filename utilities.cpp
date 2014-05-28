#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cfloat>

#include "utilities.hpp"
#include "genome.hpp"

using namespace std;

vector<Genome> obtain_genomes(string myfile_name)
{
  vector<Genome> all_Genomes;
  
  string currName, currDNA;         // Auxiliar variables
  ifstream myfile;
  myfile.open(myfile_name.c_str());
  while (myfile >> currName && myfile >> currDNA) {
    all_Genomes.push_back(Genome(currName, currDNA));
  }
  myfile.close();
  
  return all_Genomes;
}



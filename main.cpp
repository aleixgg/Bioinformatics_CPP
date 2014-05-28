#include <iostream>
#include <fstream>
#include <vector>
#include <climits>
#include <cmath>

#include "utilities.hpp"
#include "genome.hpp"
#include "pairwise_alignment.hpp"
#include "multiple_alignment.hpp"
#include "guide_tree.hpp"
 
using namespace std;

int main (int argc, char** argv)
{
  // Use original genomes downloaded from Genbank (NOT CORRECTED!)
  string genomes_filename = "genome_bad.txt";
  vector< Genome > all_genomes = obtain_genomes(genomes_filename);
  int num_genomes = all_genomes.size();

  // Create phylogenetic trees with UPGMA and Neighbor-Joining
  GuideTree upgma_tree(all_genomes), nj_tree(all_genomes);
  upgma_tree.CreateTree("upgma");
  nj_tree.CreateTree("neighbor_joining");
  
  // Output guide trees to file
  ofstream output_upgma("upgma_tree_bad.nwk"), output_nj("nj_tree_bad.nwk");
  upgma_tree.Output(output_upgma);
  nj_tree.Output(output_nj);
  output_upgma.close();
  output_nj.close();
  
  // Create multiple alignment with UPGMA guiding tree
  MultipleAlignment mult_align(upgma_tree.Output(), all_genomes);
  mult_align.Align();
  vector<Genome> solution = mult_align.Output();
  
  // Print MSA to file
  ofstream msa_output("msa_upgma_bad.txt");
  for (int i = 0; i < int(solution[0].size() / 80); ++i) {
	  for (int j = 0; j < solution.size(); ++j) {
	    msa_output << solution[j].GetName() << " ";
	    msa_output << solution[j].substr(80 * i, 80) << endl;
	  }
	  msa_output << endl;   
	}
	msa_output.close();
	
	// Use corrected data
  string filename = "genome_good.txt";
  all_genomes = obtain_genomes(filename);
  num_genomes = all_genomes.size();

  // Create phylogenetic trees with UPGMA and Neighbor-Joining
  upgma_tree = GuideTree(all_genomes);
  nj_tree    = GuideTree(all_genomes);
  upgma_tree.CreateTree("upgma");
  nj_tree.CreateTree("neighbor_joining");
  
  // Output guide trees to file
  output_upgma.open("upgma_tree_good.nwk"); 
  output_nj.open("nj_tree_good.nwk");
  upgma_tree.Output(output_upgma);
  nj_tree.Output(output_nj);
  output_upgma.close();
  output_nj.close();
  
  // Create multiple alignment 
  mult_align = MultipleAlignment(upgma_tree.Output(), all_genomes);
  mult_align.Align();
  solution = mult_align.Output();
  
  // Print MSA to file
  msa_output.open("msa_upgma_good.txt");
  for (int i = 0; i < int(solution[0].size() / 80); ++i) {
	  for (int j = 0; j < solution.size(); ++j) {
	    msa_output << solution[j].GetName() << " ";
	    msa_output << solution[j].substr(80 * i, 80) << endl;
	  }
	  msa_output << endl;   
	}
	msa_output.close();
	
  return 0;  
}

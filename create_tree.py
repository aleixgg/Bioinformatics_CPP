import sys
from Bio import Phylo

# UPGMA tree
tree_upgma = Phylo.read("upgma_tree_good.nwk", 'newick')
tree_upgma.ladderize()
Phylo.draw(tree_upgma)

tree_upgma = Phylo.read("upgma_tree_bad.nwk", 'newick')
tree_upgma.ladderize()
Phylo.draw(tree_upgma)

# NJ tree
tree_nj = Phylo.read("nj_tree_good.nwk", 'newick')
tree_nj.root_at_midpoint()
tree_nj.ladderize()
Phylo.draw(tree_nj)

tree_nj = Phylo.read("nj_tree_bad.nwk", 'newick')
tree_nj.root_at_midpoint()
tree_nj.ladderize()
Phylo.draw(tree_nj)

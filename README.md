Instructions to compile and run the code

Option 1:
  * $ make
  * $ ./Bioinformatics
  * $ python create_tree.py
  
Opció 2:
  * $ g++ -O2 main.cpp genome.cpp pairwise_alignment.cpp multiple_alignment.cpp 
              guide_tree.cpp utilities.cpp -o Bioinformatics
  * $ ./Bioinformatics
  * $ python create_tree.py
  
Important! In order to execute the python script 'create_tree.py' the Biopython
package has to be installed on the machine.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Instructions to compile and run the code

Option 1:
  * $ make
  * $ ./Bioinformatics
  * $ python create_tree.py
  
Option 2:
  * $ g++ -O2 main.cpp genome.cpp pairwise_alignment.cpp multiple_alignment.cpp 
              guide_tree.cpp utilities.cpp -o Bioinformatics
  * $ ./Bioinformatics
  * $ python create_tree.py
  
Important! In order to execute the python script 'create_tree.py' the Biopython
package has to be installed on the machine.

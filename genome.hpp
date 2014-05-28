#ifndef GENOME_HPP
#define GENOME_HPP

#include <string>

class Genome
{
public:
  Genome();
  Genome(std::string rName, std::string rDNA);
  Genome(const Genome& rGenome);
  ~Genome();
  
  Genome& operator= (const Genome&);
  char& operator[] (int i);
  int size();
  std::string substr(int i, int l);
  void resize(int l);
  
  void SetName(std::string rName);
  void SetDNA(std::string rDNA);
  std::string GetDNA();
  std::string GetName();
  
private:
  std::string mDNA;
  std::string mName;
};

#endif // GENOME_HPP

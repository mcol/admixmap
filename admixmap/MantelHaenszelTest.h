// *-*-C++-*-*
#include <vector>
#include <fstream>
class IndividualCollection;
class Genome;

class MantelHaenszelTest{
public:
  MantelHaenszelTest();
  ~MantelHaenszelTest();
  void Initialise(unsigned NumStates, unsigned NumLoci);
  void Update(const IndividualCollection* , const Genome&);
  void Output(const char* filename, const Genome& , unsigned NumIters, const std::vector<std::string>& LocusLabels);

private:
  unsigned K, Ksq;//number of states (populations or block states)
  std::vector<std::vector<unsigned> > CountTable;
  std::vector<double> Exp;
  std::vector<double> Var;

  std::ofstream outfile; 


};

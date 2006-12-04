// *-*-C++-*-*
#include "ScoreTestBase.h"
#include <vector>
#include <fstream>
class IndividualCollection;
class Genome;

class MantelHaenszelTest : public ScoreTestBase{
public:
  MantelHaenszelTest();
  ~MantelHaenszelTest();

  void Reset(){};
  void Initialise(unsigned NumStates, unsigned NumLoci);
  void Update(const IndividualCollection* , const Genome&);
  void Output(const char* filename, const Genome& , unsigned NumIters, const std::vector<std::string>& LocusLabels);

private:
  unsigned K, Ksq;//number of states (populations or block states)
  std::vector<std::vector<unsigned> > CountTable;
  std::vector<double> Score;
  std::vector<double> ScoreSq;
  std::vector<double> Info;

  std::ofstream outfile; 


};

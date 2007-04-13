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
  void Initialise(unsigned NumStates, const Genome* const Loci, const char* filename, LogWriter& Log);
  void Update(const IndividualCollection* , const Genome&);
  void Output(const char* filename,  unsigned NumIters, const std::vector<std::string>& LocusLabels, bool final);

private:
  unsigned K, Ksq;//number of states (populations or block states)
  const Genome* Loci;
  std::vector<std::vector<unsigned> > CountTable;
  std::vector<double> Score;
  std::vector<double> ScoreSq;
  std::vector<double> Info;

};

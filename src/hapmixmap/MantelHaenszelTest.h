// *-*-C++-*-*
#include "ScoreTestBase.h"
#include "bclib/RObjectWriter.h"
#include <vector>
#include <fstream>
class IndividualCollection;
class Genome;
class LogWriter;

class MantelHaenszelTest : public ScoreTestBase{
public:
  MantelHaenszelTest();
  ~MantelHaenszelTest();

  void Reset(){};
  void Initialise(unsigned NumStates, const Genome* const Loci, const std::string& ResultsDir);
  void Update(const IndividualCollection* , const Genome&);
  void Output(const std::vector<std::string>& LocusLabels);
  void WriteFinalTable(const std::string& ResultsDir, 
		       const std::vector<std::string>& LocusLabels, LogWriter& Log);

private:
  unsigned K, Ksq;//number of states (populations or block states)
  const Genome* Loci;
  std::vector<std::vector<unsigned> > CountTable;
  std::vector<double> Score;
  std::vector<double> ScoreSq;
  std::vector<double> Info;

  RObjectWriter R;
  void OutputTest( FileWriter& outfile, const std::vector<std::string>& LocusLabels, bool final);
};

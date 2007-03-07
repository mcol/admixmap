#include "HapMixIndividualTest.h"

#define REPEAT(TIMES, WHAT) for (int i = 0; i < TIMES; ++i) { WHAT }

HapMixIndividualTest::HapMixIndividualTest()
{
}

HapMixIndividualTest::~HapMixIndividualTest()
{
}

void HapMixIndividualTest::setUp()
{
  // TODO: Create a HapMixIndividual object
  
  // A mock object is created
  options = new MockOptions();
  // The mock object doesn't know that to do. We need to instruct it.
  // Adding a controller for its method named getHWTestIndicator
  MOCKPP_CONTROLLER_FOR(MockOptions, getHWTestIndicator) HWTestIndicatorCtrl (options);
  // Telling the controller which value to return.
  HWTestIndicatorCtrl.addReturnValue(true);
  // getPopulations()
  MOCKPP_CONTROLLER_FOR(MockOptions, getPopulations) getPopulationsCtrl (options);
  getPopulationsCtrl.addReturnValue(2);
  //  getgenotypesSexColumn()
  MOCKPP_CONTROLLER_FOR(MockOptions, getgenotypesSexColumn) getgenotypesSexColumnCtrl (options);
  getgenotypesSexColumnCtrl.addReturnValue(1);
  getgenotypesSexColumnCtrl.addReturnValue(1);
  getgenotypesSexColumnCtrl.addReturnValue(1);
  //  isRandomMatingModel()
  MOCKPP_CONTROLLER_FOR(MockOptions, isRandomMatingModel) isRandomMatingModelCtrl (options);
  isRandomMatingModelCtrl.addReturnValue(true);
  
  
  inputData = new MockInputData();
  // isFemale()
  MOCKPP_CONTROLLER_FOR_EXT(MockInputData, isFemale, ext) isFemaleCtrl (inputData);
  isFemaleCtrl.addReturnValue(true);
  
  vector<short unsigned int> v11;
  v11.push_back(1); v11.push_back(2);
  
  vector<vector<short unsigned int> > v12;
  v12.push_back(v11); v12.push_back(v11);
  
  vector<vector<vector<short unsigned int> > > genotypes;
  genotypes.push_back(v12);
  inputData->setGetGenotypeDefaultReturnValue(genotypes);
  
  vector<short unsigned int> v21;
  v21.push_back(2);
  inputData->setGetHapMixGenotypeDefaultReturnValues(false, v21);
  
  // A = new MockHapMixFreqs();
  haploidGenProbs = new MockFreqArray();
  diploidGenProbs = new MockFreqArray();

  genome = new MockGenome();
  MOCKPP_CONTROLLER_FOR(MockGenome, GetNumberOfChromosomes) GetNumberOfChromosomesCtrl (genome);
  GetNumberOfChromosomesCtrl.addReturnValue(1);
  MOCKPP_CONTROLLER_FOR_EXT(MockGenome, getNumberOfLoci, ext) getNumberOfLociCtrl (genome);
  REPEAT(3, getNumberOfLociCtrl.addReturnValue(1);)
  MOCKPP_CONTROLLER_FOR(MockGenome, isX_data) isX_dataCtrl (genome);
  isX_dataCtrl.addReturnValue(false);
  MOCKPP_CONTROLLER_FOR_EXT(MockGenome, GetSizeOfChromosome, ext) GetSizeOfChromosomeCtrl (genome);
  GetSizeOfChromosomeCtrl.addReturnValue(1);
  GetSizeOfChromosomeCtrl.addReturnValue(1);
  MOCKPP_CONTROLLER_FOR(MockGenome, GetLengthOfGenome) GetLengthOfGenomeCtrl (genome);
  GetLengthOfGenomeCtrl.addReturnValue(1);
  MOCKPP_CONTROLLER_FOR(MockGenome, GetNumberOfCompositeLoci) GetNumberOfCompositeLociCtrl (genome);
  REPEAT(4, GetNumberOfCompositeLociCtrl.addReturnValue(1);)
  // GetTotalNumberOfLoci()
  MOCKPP_CONTROLLER_FOR(MockGenome, GetTotalNumberOfLoci) GetTotalNumberOfLociCtrl (genome);
  GetTotalNumberOfLociCtrl.addReturnValue(1);
  // isXLocus()
  MOCKPP_CONTROLLER_FOR_EXT(MockGenome, isXLocus, ext) isXLocusCtrl (genome);
  REPEAT(1, isXLocusCtrl.addReturnValue(false);)
  
  // Above values (addReturnValue) are sufficient to create
  // a minimal Individual.
  
  genome->activate();
  options->activate();
  inputData->activate();
  haploidGenProbs->activate();
  diploidGenProbs->activate();
  
  HapMixIndividual::SetStaticMembers(
      genome, // was Loci
      options,
      *haploidGenProbs,
      *diploidGenProbs);
  
  testIndiv1 = NULL;
  try {
    testIndiv1 = new HapMixIndividual(0, options, inputData);
  } catch (string s) {
    // Catch strings that may come from an individual
    throw CppUnit::Exception(CppUnit::Message(
        "Couldn't create HapMixIndividual", "reason: " + s));
  } catch (CppUnit::Exception e) {
    // CppUnit exceptions with line numbers.
    e.setMessage(CppUnit::Message(
        "Can't create HapMixIndividual.", e.what()));
    throw e;
  }
  return;
}

void HapMixIndividualTest::tearDown()
{
  genome->verify();
  inputData->verify();
  haploidGenProbs->verify();
  diploidGenProbs->verify();
  if (testIndiv1) {
    delete testIndiv1;
  }
  delete options;
  delete inputData;
  return;
}

void HapMixIndividualTest::testUnorderedGenotypeProbs()
{
  CPPUNIT_ASSERT(testIndiv1 != NULL);
  // We need another mock genome to test this function.
  MockGenome *genome2 = new MockGenome();
  
  // Teach the mock objects what to do to allow the tested function
  // to be called.
  
  // GetNumberOfStates()
  MOCKPP_CONTROLLER_FOR_EXT(MockGenome, GetNumberOfStates, ext1) GetNumberOfStatesCtrl (genome2);
  GetNumberOfStatesCtrl.addReturnValue(2);
  // getRelativeLocusNumber()
  MOCKPP_CONTROLLER_FOR_EXT(MockGenome, getRelativeLocusNumber, ext) getRelativeLocusNumberCtrl (genome2);
  getRelativeLocusNumberCtrl.addReturnValue(1);
  // getChromosomeNumber()
  MOCKPP_CONTROLLER_FOR_EXT(MockGenome, getChromosomeNumber, constint) getChromosomeNumberCtrl (genome2);
  getChromosomeNumberCtrl.addReturnValue(false);
  // GetNumberOfCompositeLoci()
  MOCKPP_CONTROLLER_FOR(MockGenome, GetNumberOfCompositeLoci) GetNumberOfCompositeLociCtrl (genome2);
  REPEAT(0, GetNumberOfCompositeLociCtrl.addReturnValue(1); )
  // getNumberOfLoci()
  MOCKPP_CONTROLLER_FOR_EXT(MockGenome, getNumberOfLoci, ext) getNumberOfLociCtrl (genome2);
//  REPEAT(0, getNumberOfLociCtrl.addReturnValue(1);)
  
  MockChromosome *chromosome = new MockChromosome();
  // getHiddenStateProbs();
  MOCKPP_CONTROLLER_FOR_EXT(MockChromosome, getHiddenStateProbs, ext) getHiddenStateProbsCtrl (chromosome);

  vector<double> hiddenStateProbs;
  hiddenStateProbs.push_back(0.1);
  hiddenStateProbs.push_back(0.2);
  hiddenStateProbs.push_back(0.3);
  hiddenStateProbs.push_back(0.4);
  getHiddenStateProbsCtrl.addReturnValue(hiddenStateProbs);

  chromosome->activate();
  
  // Composite Locus
  MockCompositeLocus *compositeLocus;
  compositeLocus = new MockCompositeLocus();
  vector<double> condHapPairProbs;
#define HPP(V1, V2, V3, V4) \
  condHapPairProbs.clear(); \
  condHapPairProbs.push_back(V1); \
  condHapPairProbs.push_back(V2); \
  condHapPairProbs.push_back(V3); \
  condHapPairProbs.push_back(V4); \
  compositeLocus->getConditionalHapPairProbsAddReturnValue(condHapPairProbs);
  
  // Data as described in a prototype R script:
  // http://actin.ucd.ie/trac/genepi/attachment/ticket/16/nst.R
  HPP(0.1, 0.1, 0.2, 0.6)
  HPP(0.2, 0.1, 0.1, 0.6)
  HPP(0.1, 0.6, 0.2, 0.1)
  HPP(0.6, 0.1, 0.2, 0.1)
  
  genome2->setMockChromosome(chromosome);
  genome2->setMockCompositeLocus(compositeLocus);
  genome2->activate();
  
  HapMixIndividual::setGenome(genome2);
  
  int locusNo = 0;
  testIndiv1->setPopulations(2);
  
  CPPUNIT_ASSERT(Individual::getNumberOfHiddenStates() == 2);
  
  testIndiv1->calculateUnorderedGenotypeProbs(locusNo);
 
  vector<vector<double> > upr = testIndiv1->getUnorderedProbs(locusNo);
  vector<vector<double> >::iterator i;

  double sum = 0;
  // Unfortunately, we can't call `accumulate' here because of the 
  // nested vectors.
  for (i = upr.begin(); i != upr.end(); ++i) {
    sum += (*i)[0];
    CPPUNIT_ASSERT((*i)[0] >= 0);
  }
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, sum, 1e-8);

  // Check the result
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.32, upr[0][0], 1e-8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.43, upr[1][0], 1e-8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25, upr[2][0], 1e-8);
  
  genome2->verify();
  delete genome2;
  return;
}

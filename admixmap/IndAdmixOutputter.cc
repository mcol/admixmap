#include "IndAdmixOutputter.h"

using namespace std;

IndAdmixOutputter::IndAdmixOutputter(AdmixOptions* options,Genome* Loci,string* PopulationLabels)
{
  _options           = options;
  _Loci              = Loci;
  _PopulationLabels  = PopulationLabels;

  _iterations        = 0;
  _totalIndividuals  = 0;
  _currentIndividual = 0;


  _analysisType = options->getAnalysisTypeIndicator();
  _isLogistic = false;
  if((_analysisType==3)||(_analysisType==4)){
    _isLogistic = true;
  }

  _ModelIndicator = options->isRandomMatingModel();

  if (_options->getLocusForTest() >= (int)_Loci->GetNumberOfCompositeLoci()){
    cerr << "locusfortest is greater than number of loci" << endl;
    exit(0);
  }
  
  _out.open(options->getIndAdmixtureFilename(), ios::out );
  _out << "structure(.Data=c(" << endl;
}

IndAdmixOutputter::~IndAdmixOutputter()
{
  int dimOne = 0;

  if (_options->getLocusForTestIndicator()){
    dimOne += 2;
    if (_options->getPopulations() > 1 ){
      dimOne += 2;
    }
  }
  if( _analysisType > 1 ){//need this with commented bits below removed
    dimOne -= 1;
  }
  if (_options->getPopulations() > 0 && _ModelIndicator ){
    dimOne += 2 * _options->getPopulations();
  }
  else if(_options->getPopulations() > 0 && !_ModelIndicator ){
    dimOne += _options->getPopulations();
  }

  if( (_analysisType != 1) && (_analysisType != 0) ){
    dimOne++;
  }

  if( _options->getRhoIndicator() ){
     if(_options->isRandomMatingModel())
        dimOne += 2;
     else
        dimOne++;
  }
  
  _out << ")," << endl;
  _out << ".Dim = c(" << dimOne << "," << _totalIndividuals << "," << _iterations << ")," << endl;

  _out << ".Dimnames=list(c(";

  for( int i = 0; i < _options->getPopulations(); i++ ){
     _out << _PopulationLabels[i] << ",";
     if(_options->isRandomMatingModel() )
        _out << _PopulationLabels[i] << ",";
  }

  if( _options->getRhoIndicator() ){
     if(_options->isRandomMatingModel())
        _out << "\"rho0\",\"rho1\",";
     else
        _out << "\"rho\",";
  }

  if( (_analysisType < 0) ){
     _out << "\"Log-likelihood\"";
  }

//   if( _analysisType > 1 ){
//     _out << "\"ExpectedOutcomeVar\"";
//   }

  if (_options->getLocusForTestIndicator()){
    if( _options->getPopulations() > 1 ){
      _out << ",\"LocusAncestry1\",\"LocusAncestry2\"";
    }
    _out << ",\"Hap1\",\"Hap2\""; 
  }

  _out << "),character(0),character(0)))" << endl;
  _out.close();
}

void
IndAdmixOutputter::visitIndividual(Individual& ind, vector<int> _locusfortest, double LogLikelihood)
{
  for( int k = 0; k < _options->getPopulations(); k++ ){
     if(_options->isRandomMatingModel()){
       _out << ind.getAdmixtureProps()[k] << "," << ind.getAdmixtureProps()[ k + _options->getPopulations()] << ",";
    } else {
        _out << ind.getAdmixtureProps()[k] << ",";
     }
  }
  if( _options->getRhoIndicator() ){
     vector<double> rho = ind.getRho();
     if(_options->isRandomMatingModel())
        _out << rho[0] << "," << rho[1] << ",";
     else
        _out << rho[0] << ",";
  }
  
        
//   if( _analysisType > 1 ){
//     _out << expectedY << ",";
//   }

  if( (_analysisType < 0) ){
     _out << LogLikelihood << ",";
  }

  if (_options->getLocusForTestIndicator()){
    int ancestry[2];
    ind.GetLocusAncestry( _locusfortest[0], _locusfortest[1], ancestry );
     unsigned short **genotype = ind.getGenotype(_options->getLocusForTest());
     if(_options->getPopulations() > 1 ){
        _out << ancestry[0] << "," << ancestry[1] << ",";
     }
     if((*_Loci)(_options->getLocusForTest() )->GetNumberOfLoci() > 1 ){
       int hap[2] = {0,0};
       (*_Loci)(_options->getLocusForTest())->
	 SampleHapPair(hap,ind.getPossibleHapPairs(_options->getLocusForTest()),ancestry);
        if(_options->getPopulations() > 1 ){
           genotype = ind.getGenotype(_options->getLocusForTest() );
        }
        _out << hap[0] << "," << hap[1] << ",";
     } else {
        _out << genotype[0][0] << "," << genotype[0][1] << ",";
     }
  }
  _out << endl;
  _currentIndividual++;
}

void
IndAdmixOutputter::visitIndividualCollection(IndividualCollection& i)
{
  _iterations++;
  _totalIndividuals = i.getSize();
  _currentIndividual = 0;
}


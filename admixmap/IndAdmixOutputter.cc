#include "IndAdmixOutputter.h"

using namespace std;

IndAdmixOutputter::IndAdmixOutputter(const AdmixOptions* const options, const Genome* const Loci, const Vector_s& PopulationLabels)
{
  _options           = options;
  _Loci              = Loci;
  _PopulationLabels  = &PopulationLabels;

  _iterations        = 0;
  _totalIndividuals  = 0;
  _currentIndividual = 0;


  _RandomMatingModelIndicator = options->isRandomMatingModel();

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

//number of cols for admixture proportions
  if (_options->getPopulations() > 0){
    if( _RandomMatingModelIndicator ){
      dimOne += 2 * _options->getPopulations();
    }
    else {
      dimOne += _options->getPopulations();
    }
  }

  //number of cols for sumintensities
  if( !_options->isGlobalRho() ){
     if(_options->isRandomMatingModel())
        dimOne += 2;
     else
        dimOne++;
  }
//   //column for loglikelihood
//   if(_totalIndividuals == 1)
//     dimOne++;
  
  _out << ")," << endl;
  _out << ".Dim = c(" << dimOne << "," << _totalIndividuals << "," << _iterations << ")," << endl;

  _out << ".Dimnames=list(c(";

  for( int i = 0; i < _options->getPopulations(); i++ )
    {
      _out << "\""<< (*_PopulationLabels)[i];
    if(_options->isRandomMatingModel() )
      _out << "1";
    _out << "\",";
    }
  if(_options->isRandomMatingModel() )
    for( int i = 0; i < _options->getPopulations(); i++ )
      _out << "\""<<(*_PopulationLabels)[i] << "2\",";
  

  if( !_options->isGlobalRho() ){
     if(_options->isRandomMatingModel())
        _out << "\"rho1\",\"rho2\",";
     else
        _out << "\"rho\",";
  }

//   if( (_totalIndividuals == 1) ){
//      _out << "\"Log-likelihood\"";
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

void IndAdmixOutputter::visitIndividual(const AdmixedIndividual& ind, const vector<int> _locusfortest)
{
  for( int k = 0; k < _options->getPopulations(); k++ )
    _out << ind.getAdmixtureProps()[k] << ",";
  if(_options->isRandomMatingModel())
    for( int k = 0; k < _options->getPopulations(); k++ )
      _out << ind.getAdmixtureProps()[ k + _options->getPopulations()] << ",";
  
  if( !_options->isGlobalRho() ){
     vector<double> rho = ind.getRho();
     _out << rho[0] << ","; 
     if(_options->isRandomMatingModel())
        _out << rho[1] << ",";
  }
  
        
//   if( (_totalIndividuals == 1 ) ){
//      _out << LogLikelihood << ",";
//   }

  if (_options->getLocusForTestIndicator()){
    int ancestry[2];
    ind.GetLocusAncestry( _locusfortest[0], _locusfortest[1], ancestry );
    //vector<vector<unsigned short> > genotype_ = ind.getGenotype(_options->getLocusForTest());
     if(_options->getPopulations() > 1 ){
        _out << ancestry[0] << "," << ancestry[1] << ",";
     }
     //if((*_Loci)(_options->getLocusForTest() )->GetNumberOfLoci() > 1 ){
       const int* happair = ind.getSampledHapPair(_options->getLocusForTest());

       //if(_options->getPopulations() > 1 ){
       //  genotype_ = ind.getGenotype(_options->getLocusForTest() );
       //}
        _out << happair[0] << "," << happair[1] << ",";
	//} else {
        //_out << genotype_[0][0] << "," << genotype_[0][1] << ",";
	//}
  }
  _out << endl;
  _currentIndividual++;
}

void IndAdmixOutputter::visitIndividualCollection(const IndividualCollection& i)
{
  _iterations++;
  _totalIndividuals = i.getSize();
  _currentIndividual = 0;
}


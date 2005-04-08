#include "AlleleFreqOutputter.h"

using namespace std;

AlleleFreqOutputter::AlleleFreqOutputter( AdmixOptions* options, string* PopulationLabels )
{
  _options = options;
  _PopulationLabels = PopulationLabels;
  _iterations = 0;
  _numLoci = 0;

  _out.open(_options->getAlleleFreqOutputFilename(), ios::out );
  if( !_out && options->getAnalysisTypeIndicator() >= 0){
    cerr << "Warning: Couldn't open allelefreqsoutputfile: " << options->getAlleleFreqOutputFilename() << endl;
    //exit( 1 );
  }
  else
    _out << "structure(.Data=c(" << endl;
}

AlleleFreqOutputter::~AlleleFreqOutputter()
{
  _out << ")," << endl;
  _out << ".Dim = c(";
  _out << _options->getPopulations()+1 << ",";
  _out << _numLoci << ",";
  _out << _iterations;
  _out << ")," << endl;
  _out << ".Dimnames=list(c(\"Locus\",";
  for (int i=0; i<_options->getPopulations(); i++){
    _out << _PopulationLabels[i];
    if(i<_options->getPopulations()-1){
      _out << ",";
    }
  }
  _out << "), character(0), character(0)))" << endl;
  _out.close();
}

void AlleleFreqOutputter::OutputAlleleFreqs(AlleleFreqs* A)
{
  Matrix_d freqs;
  Matrix_d meanfreqs;
  for( int i = 0; i < A->GetNumberOfCompositeLoci(); i++ ){
    freqs = A->GetAlleleFreqs(i);
    // calculation of posterior means is not required
    // meanfreqs = A->GetSumAlleleFreqs(i)/_iterations;
    if(_iterations==1){
      _numLoci += freqs.GetNumberOfRows();
    }
    for( int k = 0; k < freqs.GetNumberOfRows(); k++ ){
      _out << A->getLocus(i)->GetLabel(0) << ",";
      for( int l = 0; l < _options->getPopulations(); l++ ){
	_out << freqs(k,l) << ",";
      }
      _out << endl;
    }
  }
}

void
AlleleFreqOutputter::visitChromosome(Chromosome&)
{
}

void
AlleleFreqOutputter::visitGenome(Genome&)
{
  _iterations++;
}

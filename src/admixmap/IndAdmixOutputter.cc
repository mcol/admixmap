/**
 *   ADMIXMAP
 *   IndAdmixOutputter.cc
 *   Class to output individual admix parameters
 *   Copyright (c) 2002-2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *
 * This program is free software distributed WITHOUT ANY WARRANTY.
 * You can redistribute it and/or modify it under the terms of the GNU General Public License,
 * version 2 or later, as published by the Free Software Foundation.
 * See the file COPYING for details.
 *
 */
#include "IndAdmixOutputter.h"
#include "AdmixIndividualCollection.h"


using namespace std;

using genepi::RhoType;



IndAdmixOutputter::IndAdmixOutputter(const AdmixOptions& options, const Genome& Loci, const Vector_s& PopulationLabels):
  _options(options), _Loci(Loci), _PopulationLabels(PopulationLabels),  _RandomMatingModelIndicator(options.isRandomMatingModel())
{
  _iterations        = 0;
  _totalIndividuals  = 0;
  _currentIndividual = 0;

  if (_options.getLocusForTest() >= (int)_Loci.GetNumberOfCompositeLoci()){
    cerr << "locusfortest is greater than number of loci" << endl;
    exit(0);
  }

  _out.open(options.getIndAdmixtureFilename());
}

IndAdmixOutputter::~IndAdmixOutputter(){
  //finish writing R object
  vector<vector<string> > dimnames(1);

  for( int i = 0; i < _options.getPopulations(); i++ ) {
    dimnames[0].push_back(_PopulationLabels[i]);

    if(_RandomMatingModelIndicator )
      dimnames[0][i].append( "1" );
  }

  if(_RandomMatingModelIndicator )
    for( int i = 0; i < _options.getPopulations(); i++ )
      dimnames[0].push_back( _PopulationLabels[i] + "2");


  if( !_options.isGlobalRho() ){
    if(_options.isRandomMatingModel()){
      dimnames[0].push_back("rho1");
      dimnames[0].push_back("rho2");
    }
    else
      dimnames[0].push_back("rho");
  }

  if (_options.getLocusForTestIndicator()){
    if( _options.getPopulations() > 1 ){
      dimnames[0].push_back("LocusAncestry1");
      dimnames[0].push_back("LocusAncestry2");
    }
    dimnames[0].push_back("Hap1");
    dimnames[0].push_back("Hap2");
  }

  vector<int> dims;
  dims.push_back(dimnames[0].size());
  dims.push_back(_totalIndividuals);
  dims.push_back(_iterations);
  _out.close(dims, dimnames);
}

void IndAdmixOutputter::visitIndividual(const PedBase & ind, const vector<int> _locusfortest)
{

  //output individual admixture proportions
  for( int k = 0; k < _options.getPopulations(); k++ ){
    _out << ind.getAdmixtureProps()[k] ;
  }

  //if random mating, output admixture proportions for second gamete
  if(_options.isRandomMatingModel())
    for( int k = 0; k < _options.getPopulations(); k++ )
      _out << ind.getAdmixtureProps()[ k + _options.getPopulations()] ;

  //output individual sumintensities
  if( !_options.isGlobalRho() ){
     const RhoType & rho = ind.getRho();
     //first gamete
     _out << rho[0];
     //second gamete if there is one
     if(_options.isRandomMatingModel())
       _out << rho[1] ;
  }

  if (_options.getLocusForTestIndicator()){
    int ancestry[2];
    ind.GetLocusAncestry( _locusfortest[0], _locusfortest[1], ancestry );
    //vector<vector<unsigned short> > genotype_ = ind.getGenotype(_options.getLocusForTest());
     if(_options.getPopulations() > 1 ){
        _out << ancestry[0] << ancestry[1];
     }
     //if(_Loci(_options.getLocusForTest() )->GetNumberOfLoci() > 1 ){
       const int* happair = ind.getSampledHapPair(_options.getLocusForTest());

       //if(_options.getPopulations() > 1 ){
       //  genotype_ = ind.getGenotype(_options.getLocusForTest() );
       //}
        _out << happair[0] << happair[1];
	//} else {
        //_out << genotype_[0][0] << genotype_[0][1];
	//}
  }
  _out << bclib::newline;
  _currentIndividual++;
}

void IndAdmixOutputter::visitIndividualCollection(const AdmixIndividualCollection& i)
{
  _iterations++;
  _totalIndividuals = i.getSize();
  _currentIndividual = 0;
}

/**
   \file HapMapData.cpp
   functions to write HAPMIXMAP formatted locusfiles and genotypesfiles, based on HapMap data.
   This file is part of FPHD

   This program is free software distributed WITHOUT ANY WARRANTY. 
   You can redistribute it and/or modify it under the terms of the GNU General Public License, 
   version 2 or later, as published by the Free Software Foundation. 
   See the file COPYING for details.

   Copyright (c) David O'Donnell 2007
*/
#include "HapMapData.h"
#include "HapMapLegend.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>

using namespace::std;

void OpenLocusFiles( vector<ofstream*>& locusfiles, const char* fileprefix, unsigned NumSubChrs);
void OpenGenotypeFiles(vector<ofstream*>& genofiles, const char* fileprefix, unsigned NumSubChrs);
void WriteGenotypesHeader(vector<ofstream*>& genotypesfiles, HapMapLegend& Legend, unsigned first, unsigned last);
void WriteLocusFileBody(HapMapLegend& Legend, vector<ofstream*>& locusfiles,
			unsigned first, unsigned last, bool beVerbose);
void WriteGenotypesFileBody(vector<ofstream*>& genotypesfiles, HapMapLegend& Legend, 
			    unsigned first, unsigned last, const string& prefix, bool beVerbose);
vector<unsigned> FindMonomorphicLoci(vector<unsigned> AlleleCounts[2]);
void WriteLocusFile(HapMapLegend& Legend, const char* fileprefix, 
		    unsigned first, unsigned last, bool beVerbose){

  vector<ofstream*> locusfiles;

  OpenLocusFiles(locusfiles, fileprefix, Legend.getNumSubChromosomes());
  WriteLocusFileBody(Legend, locusfiles, first, last, beVerbose);

  //clean up
  for(vector<ofstream*>::iterator i = locusfiles.begin(); i != locusfiles.end(); ++i){
    (*i)->close();
    delete *i;
  }
  locusfiles.clear();
}

void WriteGenotypesFile(HapMapLegend& Legend, const char* fileprefix, const string& dataprefix, 
			unsigned first, unsigned last, bool beVerbose){

  vector<ofstream*> genotypefiles;

  OpenGenotypeFiles(genotypefiles, fileprefix, Legend.getNumSubChromosomes());
  WriteGenotypesHeader(genotypefiles, Legend, first, last);
  WriteGenotypesFileBody(genotypefiles, Legend,  first, last, dataprefix, beVerbose);

  //clean up
  for(vector<ofstream*>::iterator i = genotypefiles.begin(); i != genotypefiles.end(); ++i){
    (*i)->close();
    delete *i;
  }
  genotypefiles.clear();
  //  genotypefile.close();
}

void OpenLocusFiles(vector<ofstream*>& locusfiles, const char* fileprefix, unsigned NumSubChrs){
  for(unsigned j = 0; j < NumSubChrs; ++j){
    stringstream locusfilename;
    if(NumSubChrs==1)
      locusfilename << fileprefix << ".txt";
    else
      locusfilename << fileprefix << "_" << j+1 << ".txt";
    
    locusfiles.push_back(new ofstream(locusfilename.str().c_str()));
    if (!locusfiles[j]->is_open()) {
      cout << "Error: cannot open locusfile\n";
      exit(1);
    }
    *(locusfiles[j]) << "\"SNPid\"\t\"NumAlleles\"\t\"DistanceinMb\"\n"; //locusfile header
    *(locusfiles[j]) << setiosflags(ios::fixed) << setprecision(8);
  }
}

void WriteLocusFileBody(HapMapLegend& Legend, vector<ofstream*>& locusfiles, 
			unsigned first, unsigned last, bool beVerbose){

  unsigned long prev = 0, position = 0;
  vector<bool>isFirstLocus (locusfiles.size(), true);
  for(unsigned locus = first; locus < last; ++locus){
    if (beVerbose)
      cout << "\r" << locus + 1 << " loci" << flush;

    prev = position;
    position = Legend[locus].position;
    string SNPID = Legend.getRSNumber(locus);
    if ((position - prev) > 0) {    //strictly greater than to avoid having comp loci
      for(vector<unsigned>::const_iterator j = Legend[locus].subchr.begin(); 
	  j != Legend[locus].subchr.end(); ++j){
	*(locusfiles[*j]) << SNPID << "\t" << 2 << "\t";
	if (isFirstLocus[*j])
	  *(locusfiles[*j]) << "#"; //missing value for first distance on chr
	else
	  //converts position in basepairs to distance in centiMorgans 
	  //locusfile << 0.0000013 * (position-prev);
	  //convert position to distance in megabases
	  *(locusfiles[*j]) << 0.000001 * (position - prev);
	*(locusfiles[*j]) << /*chrnumber << */ endl;
	isFirstLocus[*j] = false;
      }

    }//end if position-prev >0
  }
  if (beVerbose)
    cout << endl;
}

void OpenGenotypeFiles(vector<ofstream*>& genotypesfiles, const char* fileprefix, unsigned NumSubChrs){
  for(unsigned j = 0; j < NumSubChrs; ++j){
    stringstream genofilename;
    if(NumSubChrs==1)
      genofilename << fileprefix << ".txt";
    else
      genofilename << fileprefix << "_" << j+1 << ".txt";
    
    genotypesfiles.push_back(new ofstream(genofilename.str().c_str()));
    if (!genotypesfiles[j]->is_open()) {
      cout << "Error: cannot open genotypesfile\n";
      exit(1);
    }
    //Write genotypesfile header  
    *(genotypesfiles[j]) << "\"Gameteid\"";
    //*(genofiles[j]) << setiosflags(ios::fixed) << setprecision(8);
  }
}

void WriteGenotypesHeader(vector<ofstream*>& genotypesfiles, HapMapLegend& Legend, unsigned first, unsigned last){
  for(unsigned locus = first; locus < last; ++locus)
    for(vector<unsigned>::const_iterator j = Legend[locus].subchr.begin(); 
	j != Legend[locus].subchr.end(); ++j){
      *(genotypesfiles[*j]) << "\t" << Legend.getRSNumber(locus);
    }

  for(unsigned j = 0; j < Legend.getNumSubChromosomes(); ++j)
    *(genotypesfiles[j]) << endl;
}

void WriteGenotypesFileBody(vector<ofstream*>& genotypesfiles, HapMapLegend& Legend, 
			    unsigned first, unsigned last, const string& prefix, bool beVerbose){
  //Note: sample file is arranged with founders first, then children.
  //children are omitted from phased file
  string ID, scrap;

  //open phased file to use to loop through gametes
  ifstream phasedfile((prefix + "_phased.txt").c_str());
  if (!phasedfile.is_open()) {
    cout << "Could not open phased file " << prefix << "_phased.txt" << endl;
    exit(1);
  }
  ifstream samplefile((prefix + "_sample.txt").c_str());
  if (!samplefile.is_open()) {
    cout << "Could not open sample file " << prefix << "_sample.txt" << endl;
    exit(1);
  }

  int gamete = 0;
  int allele = 0;
  const unsigned NumLoci = Legend.size();
  //vector to counts alleles in order to check for monomorphic loci
  vector<unsigned> NullCounts(NumLoci, 0);
  vector<unsigned> AlleleCounts[2] = {NullCounts, NullCounts};
  phasedfile >> allele;
  while (!phasedfile.eof()) {
    //if (beVerbose)
    //cout << "\nChromosome " << options.getChrNum() << "  " << flush;
    string suffix;//to give gametes unique IDs
    if (!(gamete % 2)){
      samplefile >> ID >> scrap;
      suffix = "_1";
    }
    else
      suffix = "_2";
    if (beVerbose)
      cout << "\r" << gamete + 1 << " " << " gametes" << " " << flush;

    //write indiv id and sex
    for(unsigned j = 0; j < Legend.getNumSubChromosomes(); ++j)
      *(genotypesfiles[j]) << ID  <<  suffix << "\t";        // << sex[indiv] <<"\t";

    //write genotypes for first chromosome
    for (unsigned j = 0; j <  NumLoci; ++j) {

      //if (j < options.getMaxLoci())
      if(j >= first && j<last){
	for(vector<unsigned>::const_iterator sub = Legend[j].subchr.begin(); 
	    sub != Legend[j].subchr.end(); ++sub){
	  *(genotypesfiles[*sub]) << allele + 1 << " ";
	}
	AlleleCounts[allele][j-first]++;
      }
      phasedfile >> allele;
    }
    ++gamete;
    for(unsigned j = 0; j < Legend.getNumSubChromosomes(); ++j)
      *(genotypesfiles[j]) << endl;
  }
  phasedfile.close();
  samplefile.close();

  vector<unsigned> MMLoci = FindMonomorphicLoci(AlleleCounts);
  if(MMLoci.size()){
    //if(beVerbose)
      cout << "\nWARNING: " << MMLoci.size() << " monomorphic loci";
  }

}

vector<unsigned> FindMonomorphicLoci(vector<unsigned> AlleleCounts[2]){
  vector<unsigned> MonomorphicIndices;
  for(unsigned locus = 0; locus < AlleleCounts[0].size(); ++locus){
    if(AlleleCounts[0][locus] == 0 || AlleleCounts[1][locus] == 0)
      MonomorphicIndices.push_back(locus);
  }
  return MonomorphicIndices;
}

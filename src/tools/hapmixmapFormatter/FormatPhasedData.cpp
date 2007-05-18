/**
 * FormatPhasedData.cpp
 * main source file for FPHD (Format Phased HapMap Data) 
 *
 * This program converts  **phased** hapmap data to HAPMIXMAP format.
 * Optionally formats a case-control genotypes file 
 *
 * NB: input data file must be named:
 * "chr#_phased.txt",
 * "chr#_sample.txt" and
 * "chr#_legend.txt",
 *
 * where # is a number from 1 to 22, and all be located in the same directory
 *
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 *
 * Copyright (c) David O'Donnell 2006, 2007
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include "FPHDOptions.h"
#include "GenotypeEncoding.h"

unsigned NUMALLELES = 2;

using namespace::std;


void WriteLocusFile(HapMapLegend& Legend, ofstream& locusfile, unsigned first, unsigned last, bool beVerbose){

  unsigned long prev = 0, position = 0;
  for(unsigned locus = first; locus <= last; ++locus){
    if (beVerbose)
      cout << "\r" << locus + 1 << " loci" << flush;

    prev = position;
    position = Legend[locus].position;
    string SNPID = Legend.getRSNumber(locus);
    if ((position - prev) > 0) {    //strictly greater than to avoid having comp loci
      locusfile << SNPID << "\t" << 2 << "\t";
      if (locus == first)
        locusfile << "#"; //missing value for first distance on chr
      else
        //converts position in basepairs to distance in centiMorgans 
        //locusfile << 0.0000013 * (position-prev);
        //convert position to distance in megabases
        locusfile << 0.000001 * (position - prev);
      locusfile << /*chrnumber << */ endl;

    }//end if position-prev >0
  }
  if (beVerbose)
    cout << endl;
}

int main(int argc, char **argv)
{
//   if (argc < 3) {
//     FPHDOptions::PrintHelpText();
//     exit(0);
//   }

  ofstream genotypesfile;
  ofstream locusfile;
  FPHDOptions options(argc, argv, genotypesfile, locusfile);
  unsigned first;//first locus to print on each chr
  unsigned last; //last    "      "      "     "
  double position = 0.0, prev = 0.0;
  string scrap;

  unsigned long NUMLOCI;
  unsigned long TOTALLOCI = 0;

  if(options.Verbose()){
    cout << "******************************" << endl
	 << "Beginning formatting" << endl;
  }

  // *** PHASE 1: write locus file and header of genotypesfile, count number of loci  ***
  genotypesfile << "\"Gameteid\"";
  locusfile << "\"SNPid\"\t\"NumAlleles\"\t\"DistanceinMb\"\n"; //locusfile header
  locusfile << setiosflags(ios::fixed) << setprecision(8);

  ifstream legendfile;
  if(options.Verbose())
    cout << "\nprocessing chromosome " << options.getChrNum() << "  " << endl;

  string prefix;
  {
    stringstream ss;
    ss << options.getPrefix() << "/chr" << options.getChrNum();
    prefix = ss.str();
  }

  //open HapMap legend file
  legendfile.clear();         //to clear fail status at eof
  legendfile.open((prefix + "_legend.txt").c_str());
  if (!legendfile.is_open()) {
    cout << "Could not open legend file " << prefix << "_legend.txt" << endl;
    exit(1);
      }
  HapMapLegend Legend(legendfile);
  
  if(options.WriteCCFile()){
    if(options.Verbose())
      cout << "Writing case-control genoypes to " << options.getOutCCFilename() << endl;
    unsigned NumTypedInd = EncodeGenotypes(Legend, options.getInCCFilename(), options.getOutCCFilename(), options.getMissingChar());
    if(options.Verbose())
      cout << NumTypedInd << " case-control individuals" << endl
	   << "first = " << Legend.getFirst() << "(" << Legend.getFirstIndex() << "), last = "
	   << Legend.getLast() << "(" << Legend.getLastIndex() << ")" << endl;
    
    Legend.OffsetLimits(options.getFlankLength());
    if(options.Verbose())
      cout << "first = " << Legend.getFirst() << "(" << Legend.getFirstIndex() << "), last = "
	   << Legend.getLast() << "(" << Legend.getLastIndex() << ")" << endl;
  }
  
  if(!options.WriteCCFile() && options.LimitedLoci()){
    first = 0;
    last = options.getNumUserLoci();
  }
  else{
    first = Legend.getFirstIndex();
    last = Legend.getLastIndex();
  }
  WriteLocusFile(Legend, locusfile, first, last, options.Verbose());
  legendfile.close();
  
  //Write genotypesfile header
  for(unsigned i = first; i <= last; ++i)
    genotypesfile << "\t" << Legend.getRSNumber(i);
  
  NUMLOCI = Legend.size();   //count number of loci
  TOTALLOCI += last - first +1;//Legend.size();
  Legend.clear();
  
  locusfile.close();
  if(options.Verbose())
    cout << "Finished writing locusfile. " << TOTALLOCI << " loci" << endl << endl;
  
  genotypesfile << endl;

  // ************************** PHASE 2: Write genotypesfile  *********************
  //Note: sample file is arranged with founders first, then children.
  //children are omitted from phased file
  string ID;

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
  phasedfile >> allele;
  while (!phasedfile.eof()) {
    //if (options.Verbose())
    //cout << "\nChromosome " << options.getChrNum() << "  " << flush;
    string suffix;//to give gametes unique IDs
    if (!(gamete % 2)){
      samplefile >> ID >> scrap;
      suffix = "_1";
    }
    else
      suffix = "_2";
    if (options.Verbose())
      cout << "\r" << gamete + 1 << " " << " gametes" << " " << flush;
    //write indiv id and sex
    genotypesfile << ID  + suffix << "\t";        // << sex[indiv] <<"\t";
    //write genotypes for first chromosome
    for (unsigned j = 0; j < NUMLOCI; ++j) {

      //if (j < options.getNumUserLoci())
      if(j >= first && j<=last)
        genotypesfile << allele + 1 << " ";
      phasedfile >> allele;
    }
    ++gamete;
    genotypesfile << endl;
  }
  phasedfile.close();
  genotypesfile.close();
  //  outcomefile.close();
  if(options.Verbose()){
    cout << endl << "Finished writing genotypesfile" << endl;
    //cout << gamete << " gametes" << endl;
  }

  if(options.Verbose()){
    cout << "Formatting complete" << endl
	 << "******************************" << endl;
  }
  return 0;
}

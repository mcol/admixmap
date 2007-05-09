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
      cout << "\rLocus    " << locus + 1 << flush;

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
}

int main(int argc, char **argv)
{
  if (argc < 3) {
    FPHDOptions::PrintHelpText();
    exit(0);
  }

  ofstream genotypesfile;
  ofstream locusfile;
  FPHDOptions options(argc, argv, genotypesfile, locusfile);
  vector<unsigned> first;//first locus to print on each chr
  vector<unsigned> last; //last    "      "      "     "
  double position = 0.0, prev = 0.0;
  string scrap;

  vector < unsigned long >NUMLOCI;
  unsigned long TOTALLOCI = 0;

  // ************************** PHASE 1: write locus file and header of genotypesfile, count number of loci  *********************
  genotypesfile << "\"Gameteid\"";
  locusfile << "\"SNPid\"\t\"NumAlleles\"\t\"DistanceinMb\"\n"; //locusfile header
  locusfile << setiosflags(ios::fixed) << setprecision(8);

  ifstream legendfile;
  for (unsigned chr = options.getFirstChr(); chr <= options.getLastChr(); ++chr) {
    cout << "\nprocessing chromosome " << chr << "  " << endl;
    stringstream ss;

    /// unzip data
      //ss << "gzip -d chr" << chr << ".gz";
      //cout << "\nUnzipping..."<<flush;
      //system(ss.str().c_str());
      //ss.clear();
      ///open HapMap legend file
      /// NB: files should be named chr1_legend, chr2_legend etc
      ss << options.getPrefix() << "/chr" << chr << "_legend.txt";
      legendfile.clear();         //to clear fail status at eof
      legendfile.open(ss.str().c_str());
      if (!legendfile.is_open()) {
        cout << "Could not open legend file " << ss.str() << endl;
        exit(1);
      }
      HapMapLegend Legend(legendfile);

      if(options.WriteCCFile()){
        EncodeGenotypes(Legend, options.getInCCFilename(), options.getOutCCFilename());
        if(options.Verbose())
          cout << "first = " << Legend.getFirst() << "(" << Legend.getFirstIndex() << "), last = "
               << Legend.getLast() << "(" << Legend.getLastIndex() << ")" << endl;

        Legend.OffsetLimits(options.getFlankLength());
        if(options.Verbose())
          cout << "first = " << Legend.getFirst() << "(" << Legend.getFirstIndex() << "), last = "
               << Legend.getLast() << "(" << Legend.getLastIndex() << ")" << endl;
      }
      //TODO:message saying outputting to this file

      if(!options.WriteCCFile() && options.LimitedLoci()){
        first.push_back(0);
        last.push_back(options.getNumUserLoci());
      }
      else{
        first.push_back(Legend.getFirstIndex());
        last.push_back(Legend.getLastIndex());
      }
      WriteLocusFile(Legend, locusfile, first[chr - options.getFirstChr()], last[chr - options.getFirstChr()], options.Verbose());
      legendfile.close();

      //Write genotypesfile header
      for(unsigned i = first[chr - options.getFirstChr()]; i <= last[chr - options.getFirstChr()];++i)
        genotypesfile << "\t" << Legend.getRSNumber(i);

      cout << endl << last[chr - options.getFirstChr()] - first[chr - options.getFirstChr()] +1 << " Loci ";
      NUMLOCI.push_back(Legend.size());   //count number of loci
      TOTALLOCI += last[chr - options.getFirstChr()] - first[chr - options.getFirstChr()] +1;//Legend.size();
      Legend.clear();
  }                             //end chr loop
  locusfile.close();
  cout << "\nFinished writing locusfile. " << TOTALLOCI << " loci" << endl;
  cout << endl;
  //}

  genotypesfile << endl;

  // ************************** PHASE 2: Write genotypesfile  *********************
  //TODO; tidy this section by moving some code into functions

  //Note: sample file is arranged with founders first, then children
  //children are omitted from phased file
  string ID;
  //       if(indiv%2) outcomefile << "#\n";
  //       else outcomefile << 9 << endl;

  unsigned abslocus = 0;

  //open first phased file to use to loop through gametes
  //NOTE: assuming the same set of individuals for all chromosomes
  stringstream ss;
  ss << options.getPrefix() << "/chr" << options.getFirstChr() << "_phased.txt";
  ifstream phasedfile(ss.str().c_str());
  if (!phasedfile.is_open()) {
    cout << "Could not open phased file " << ss.str() << endl;
    exit(1);
  }
  ss.clear();
  ss.str("");
  ss << options.getPrefix() << "/chr" << options.getFirstChr() << "_sample.txt";
  ifstream samplefile(ss.str().c_str());
  if (!samplefile.is_open()) {
    cout << "Could not open sample file " << ss.str() << endl;
    exit(1);
  }

  int gamete = 0;
  int allele = 0;
  ifstream phasedfile2;
  phasedfile >> allele;
  while (!phasedfile.eof()) {
    if (options.Verbose())
      cout << "\nChromosome " << options.getFirstChr() << "  " << flush;
    string suffix;
    if (!(gamete % 2)){
      samplefile >> ID >> scrap;
      suffix = "_1";
    }
    else
      suffix = "_2";
    if (options.Verbose())
      cout << "\n" << gamete + 1 << " " << ID + suffix << " " << flush;
    //write indiv id and sex
    genotypesfile << ID  + suffix << "\t";        // << sex[indiv] <<"\t";
    //write genotypes for first chromosome
    for (unsigned j = 0; j < NUMLOCI[0]; ++j) {

      //if (j < options.getNumUserLoci())
      if(j >= first[0] && j<=last[0])
        genotypesfile << allele + 1 << " ";
      phasedfile >> allele;
    }
    //now loop through other chromosomes 
    for (unsigned chr = options.getFirstChr() + 1; chr <= options.getLastChr(); ++chr) {
      if (options.Verbose())
        cout << "\nChromosome " << chr << "  " << flush;

      //        ss.clear();
      //        ss << options.getPrefix()<< "/chr" << chr << "_sample.txt";
      //        ifstream samplefile(ss.str().c_str());
      //        if(!samplefile.is.open()){
      //         cout << "Could not open sample file " << ss.str() <<endl;
      //         exit(1); 
      //        }

      //open next phased file
      ss.clear();
      ss << options.getPrefix() << "/chr" << chr << "_phased.txt";
      phasedfile2.clear();
      phasedfile2.open(ss.str().c_str());
      if (!phasedfile2.is_open()) {
        cout << "Could not open phased file " << ss.str() << endl;
        exit(1);
      }
      //write genotypes for this chromosome
      for (unsigned j = 0; j < NUMLOCI[chr - options.getFirstChr()]; ++j) {
        phasedfile2 >> allele;
        //if (j < options.getNumUserLoci())
        if(j >= first[chr - options.getFirstChr()] && j <= last[chr - options.getFirstChr()])
          genotypesfile << allele + 1;
      }
      phasedfile2.close();
    }
    ++gamete;
    genotypesfile << endl;
  }
  phasedfile.close();
  genotypesfile.close();
  //  outcomefile.close();
  cout << endl << "Finished writing genotypesfile" << endl;
  cout << gamete << " gametes" << endl;
}

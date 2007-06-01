/**
 * \file FormatPhasedData.cpp
 *  main source file for FPHD (Format Phased HapMap Data) .
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
#include <vector>
#include <string>
#include <sstream>
#include <iterator>
#include "FPHDOptions.h"
#include "HapMapLegend.h"
#include "HapMapData.h"
#include "UserGenotypes.h"
#include "InitialValues.h"

///number of alleles (fixed at 2)
unsigned NUMALLELES = 2;

using namespace::std;

int main(int argc, char **argv){  
  FPHDOptions options(argc, argv);
  unsigned first;//first locus to print on each chr
  unsigned last; //last    "      "      "     "

  if(options.Verbose()){
    cout << "******************************" << endl
	 << "Beginning formatting" << endl;
    cout << "\nprocessing chromosome " << options.getChrNum() << "  " << endl;
  }

  try{
    string prefix;
    {
      stringstream ss;
      ss << options.getPrefix() << "/chr" << options.getChrNum();
      prefix = ss.str();
    }
    
    // *** PHASE 1: read legend file  ***
    
    //read HapMap legend file
    HapMapLegend Legend((prefix + "_legend.txt").c_str());
    //set upper limit on loci to the total number of loci, if none specified
    options.setMaxLoci(Legend.getLastIndex());
    
    // *** PHASE 2: read case-control genotypes and set chromosome bounds  ***  
    UserGenotypes UG(Legend, options.getInCCFilename());
    if(options.WriteCCFile()){

      //set first and last loci in Legend
      Legend.setLimits(UG.getFirstTypedLocus(), UG.getLastTypedLocus());
      // 	   << "first = " << Legend.getFirst() << "(" << Legend.getFirstIndex() << "), last = "
      // 	   << Legend.getLast() << "(" << Legend.getLastIndex() << ")" << endl;
      
      //set flanking region
      Legend.OffsetLimits(options.getFlankLength());
      options.setMaxLoci(Legend.getLastIndex() - Legend.getFirstIndex());
      //     if(options.Verbose())
      //       cout << "first = " << Legend.getFirst() << "(" << Legend.getFirstIndex() << "), last = "
      // 	   << Legend.getLast() << "(" << Legend.getLastIndex() << ")" << endl;
    }
    //divide up chromosome as required
    Legend.DetermineCutPoints(options.getMaxLoci(), options.getMinOverlap());
    
    if(options.Verbose() && Legend.getNumSubChromosomes() > 1){
      cout << Legend.getNumSubChromosomes() << " sub-chromosomes" << endl;
      cout << "with sizes: ";
      for(vector<unsigned>::const_iterator n = Legend.getSubChromosomeSizes().begin(); n != Legend.getSubChromosomeSizes().end(); ++n)
	cout << *n << " ";
      cout << endl;
    }

    if(options.WriteCCFile()){
      if(Legend.getNumSubChromosomes() > 1)
	Legend.AdjustFlankingRegions(UG.getTypedLoci(), options.getFlankLength());
      if(options.Verbose())
	cout << "Writing case-control genoypes to " << options.getOutCCFilename() << endl;
      //format case-control file
      unsigned NumTypedInd = UG.FormatUserGenotypes(Legend,  
						options.getOutCCFilename(), options.getMissingChar());
      if(options.Verbose()){
	cout << NumTypedInd << " case-control individuals" << endl
	     << UG.getNumberOfTypedLoci() << " typed loci" << endl;
	//TODO: give numbers of typed loci on each sub-chromosome
      }
    }
    
    //   if(!options.WriteCCFile() /*&& options.LimitedLoci()*/){
    //     first = 0;
    //     last = options.getMaxLoci();
    //   }
    //   else{
    first = Legend.getFirstIndex();
    last = Legend.getLastIndex();
    //  }
    
    // *** PHASE 3: write locus file, count number of loci  ***
    WriteLocusFile(Legend, options.getLocusFilename(), first, last, options.Verbose());
    if(options.Verbose())
      cout << "Finished writing locusfile. "// << last - first << " loci" 
	   << endl << endl;
    
    // *** PHASE 4: Write genotypesfile, count gametes  ********
    WriteGenotypesFile(Legend, options.getGenotypesFilename(), prefix, first, last, options.Verbose());
    if(options.Verbose()){
      cout << endl << "Finished writing genotypesfile" << endl;
      //cout << gamete << " gametes" << endl;
    }
    
    // *** PHASE 5: Split initial value files  *********
    SplitInitialArrivalRateFile(options.getInitialArrivalRateFilename(), Legend);
    SplitInitialMixturePropsFile(options.getInitialMixturePropsFilename(), Legend);
    SplitInitialAlleleFreqFile(options.getInitialAlleleFreqFilename(), Legend);
    SplitInitialFreqPriorFile(options.getInitialFreqPriorFilename(), Legend);

    Legend.clear();
    
  }catch(string s){
    cerr << s;
    exit(1);
  }
  catch(...){
    cerr << "Unknown exception occurred. Contact authors for assistance\n";
    exit(1);
  }
  if(options.Verbose()){
    cout << "Formatting complete" << endl
	 << "******************************" << endl;
  }
  return 0;
}

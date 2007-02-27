/**
 * prog to convert  **phased** hapmap data to HAPMIXMAP format
 * supply chromosome number as arg. 0 (default) means all.  also requires 
 * data directory (Eur, Afr or Asian).  Optionally supply genotypes file 
 * and locus file names (default to genotypes.txt) loci.txt
 * NB: input data file must be named:
 * "chr#_phased.txt",
 * "chr#_sample.txt" and
 * "chr#_legend.txt",
 *
 * where # is a number from 1 to 22, and all be located in the data directory
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
#include <getopt.h>

#include "Formatting.h"

#define MAXCHROMOSOMES 22//maximum number of chromosomes
unsigned NUMALLELES = 2;
bool beVerbose = false;

using namespace::std;

void PrintHelpText(char **argv)
{
  cout << "----------------------------------------------------------------------------"
       << endl
       << "This program converts phased HapMap data to HAPMIXMAP format" << endl
       << "Copyright (c) David O'Donnell 2007" << endl
       << "All parts of this program are freely distributable" << endl << endl
       << "Usage: " << argv[0] << " -c<chr> [options]" << endl << endl
       << "Options: " << endl
       << "-c<chr>            - chromosome number. -c0 converts all chromosomes." << endl
       << "-v                - be verbose" << endl
       << "-p <prefix>        - prefix where HapMap files are located, defaults to '.'" << endl 
       << "-g <genotypesfile> - output genotypes file, defaults to 'genotypes.txt'" << endl
       << "-l <locusfile>     - output locus file, defaults to loci.txt" << endl
       << "-n <numloci>       - maximum number of loci to use per chromosome. Ignored with -i" << endl
       << "-i <input case-control file> - raw case-control genotypes file. "
       << "File will be formatted and used to restrict the range of HapMap loci output." << endl
       << "-o <output case-control-file>  - valid only with -i, defaults to 'CaseControlGenotypes.txt'." << endl
       << "Note that HapMap input files must be named 'chrN_phased.txt', 'chrN_sample.txt'" << endl
       << "and 'chrN_legend.txt', where N is an integer from 1 to 22"
       << endl << endl;
}

void WriteLocusFile(HapMapLegend& Legend, ofstream& locusfile, unsigned first, unsigned last){

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
    PrintHelpText(argv);
    exit(0);
  }
  int ich;
  unsigned FirstChr = 0;          //chromosome number
  string popprefix = ".";
  ofstream genotypesfile;
  ofstream locusfile;
  char* incasecontrolfilename = 0;
  char* outcasecontrolfilename = 0;
  unsigned long userloci = 1000000000;  //some number > number of HapMap loci
  bool LimitLoci = false;
  vector<unsigned> first;//first locus to print on each chr
  vector<unsigned> last; //last    "      "      "     "

  while ((ich = getopt(argc, argv, "hc:g:l:p:i:o:n:qv")) != EOF) {
    switch (ich) {
    case 'h':{
      PrintHelpText(argv);
    }
    case 'c':{
      FirstChr = atoi(optarg);
      break;
    }
    case 'p':{
      popprefix = optarg;
      break;
    }
    case 'g':{
      genotypesfile.open(optarg);
      if (!genotypesfile.is_open()) {
        cout << "Error: cannot open genotypesfile\n";
        exit(1);
      }
      cout << "Writing genotypes to " << optarg << endl;
      break;
    }
    case 'l':{
      locusfile.open(optarg);
      if (!locusfile.is_open()) {
        cout << "Error: cannot open locusfile\n";
        exit(1);
      }
      cout << "Writing locusfile to " << optarg << endl;
      break;
    }
    case 'i':{
      incasecontrolfilename = optarg;
      break;
    }
    case 'o':{
      outcasecontrolfilename = optarg;
      break;
    }

    case 'n':{                 //number of loci
      const unsigned temp = atoi(optarg);
      if (temp > 0){
        userloci = temp - 1;
        LimitLoci = true;
      }
      break;
    }
    case 'v':{                 //verbose mode
      beVerbose = true;
      break;
    }
    default:{
      cout << "Invalid args specified: " << ich << endl;
      exit(1);
    }

    }
  }
  if (!genotypesfile.is_open())
    genotypesfile.open("genotypes.txt");
  if (!locusfile.is_open())
    locusfile.open("loci.txt");
  //set default ouput ccgenotypesfilename if none specified
  if(incasecontrolfilename && !outcasecontrolfilename)
    outcasecontrolfilename = "CaseControlGenotypes.txt";
  //ofstream outcomefile("outcome_halfmissing.txt");

  double position = 0.0, prev = 0.0;
  string scrap;
  // unsigned chromosome = 0;
  unsigned lastchr = MAXCHROMOSOMES;
  if (FirstChr == 0)//all chromosomes (-c0 option)
    FirstChr = 1;
  else//only one chromosome
    lastchr = FirstChr;
  vector < unsigned long >NUMLOCI;
  unsigned long TOTALLOCI = 0;

  // ************************** PHASE 1: write locus file and header of genotypesfile, count number of loci  *********************
  genotypesfile << "\"Gameteid\"";
  locusfile << "\"SNPid\"\t\"NumAlleles\"\t\"DistanceinMb\"\n"; //locusfile header
  locusfile << setiosflags(ios::fixed) << setprecision(8);

  ifstream legendfile;
  for (unsigned chr = FirstChr; chr <= lastchr; ++chr) {
    cout << "\nChromosome " << chr << "  " << endl;
    stringstream ss;

    /// unzip data
      //ss << "gzip -d chr" << chr << ".gz";
      //cout << "\nUnzipping..."<<flush;
      //system(ss.str().c_str());
      //ss.clear();
      ///open HapMap legend file
      /// NB: files should be named chr1_legend, chr2_legend etc
      ss << popprefix << "/chr" << chr << "_legend.txt";
      legendfile.clear();         //to clear fail status at eof
      legendfile.open(ss.str().c_str());
      if (!legendfile.is_open()) {
        cout << "Could not open legend file " << ss.str() << endl;
        exit(1);
      }
      //map<string, LocusInfo> Legend;
      //map<unsigned, string> LegendIndex;
      //vector<LocusInfo> Legend;
      //ReadLegend(Legend, legendfile);
      HapMapLegend Legend(legendfile);

      //TODO: call EncodeGenotypes here
      if(incasecontrolfilename){
        EncodeGenotypes(Legend, incasecontrolfilename, outcasecontrolfilename);
        cout << "first = " << Legend.getFirst() << "(" << Legend.getFirstIndex() << "), last = "
             << Legend.getLast() << "(" << Legend.getLastIndex() << ")" << endl;

        Legend.OffsetLimits();
        cout << "first = " << Legend.getFirst() << "(" << Legend.getFirstIndex() << "), last = "
             << Legend.getLast() << "(" << Legend.getLastIndex() << ")" << endl;
      }
      //TODO:message saying outputting to this file

      //Write genotypesfile header
      //for(map<unsigned, string>::const_iterator p = LegendIndex.begin(); p != LegendIndex.end(); ++p){
      //for(vector<LocusInfo>::const_iterator p = Legend.begin(); p !=Legend.end(); ++p){
      // genotypesfile << "\t" << p->rsnumber;
      //}
      for(unsigned i = 0; i < Legend.size();++i)
        genotypesfile << "\t" << Legend.getRSNumber(i);

      if(!incasecontrolfilename && LimitLoci){
        first.push_back(0);
        last.push_back(userloci);
      }
      else{
        first.push_back(Legend.getFirstIndex());
        last.push_back(Legend.getLastIndex());
      }
      WriteLocusFile(Legend, locusfile, first[chr - FirstChr], last[chr - FirstChr]);

      legendfile.close();
      cout << endl << last[chr - FirstChr] - first[chr - FirstChr] +1 << " Loci ";
      NUMLOCI.push_back(Legend.size());   //count number of loci
      TOTALLOCI += last[chr - FirstChr] - first[chr - FirstChr] +1;//Legend.size();
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
  ss << popprefix << "/chr" << FirstChr << "_phased.txt";
  ifstream phasedfile(ss.str().c_str());
  if (!phasedfile.is_open()) {
    cout << "Could not open phased file " << ss.str() << endl;
    exit(1);
  }
  ss.clear();
  ss.str("");
  ss << popprefix << "/chr" << FirstChr << "_sample.txt";
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
    if (beVerbose)
      cout << "\nChromosome " << FirstChr << "  " << flush;
    if (!(gamete % 2))
      samplefile >> ID >> scrap;
    else
      ID.append("_2");
    if (beVerbose)
      cout << "\n" << gamete + 1 << " " << ID << " " << flush;
    //write indiv id and sex
    genotypesfile << ID << "\t";        // << sex[indiv] <<"\t";
    //write genotypes for first chromosome
    for (unsigned j = 0; j < NUMLOCI[0]; ++j) {

      //if (j < userloci)
      if(j >= first[0] && j<=last[0])
        genotypesfile << allele + 1 << " ";
      phasedfile >> allele;
    }
    //now loop through other chromosomes 
    for (unsigned chr = FirstChr + 1; chr <= lastchr; ++chr) {
      if (beVerbose)
        cout << "\nChromosome " << chr << "  " << flush;

      //        ss.clear();
      //        ss << popprefix<< "/chr" << chr << "_sample.txt";
      //        ifstream samplefile(ss.str().c_str());
      //        if(!samplefile.is.open()){
      //         cout << "Could not open sample file " << ss.str() <<endl;
      //         exit(1); 
      //        }

      //open next phased file
      ss.clear();
      ss << popprefix << "/chr" << chr << "_phased.txt";
      phasedfile2.clear();
      phasedfile2.open(ss.str().c_str());
      if (!phasedfile2.is_open()) {
        cout << "Could not open phased file " << ss.str() << endl;
        exit(1);
      }
      //write genotypes for this chromosome
      for (unsigned j = 0; j < NUMLOCI[chr - FirstChr]; ++j) {
        phasedfile2 >> allele;
        //if (j < userloci)
        if(j >= first[chr - FirstChr] && j <= last[chr - FirstChr])
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

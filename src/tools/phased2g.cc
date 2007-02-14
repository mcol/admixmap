/**
 * prog to convert  **phased** hapmap data to ADMIXMAP format
 *
 * Copyright (c) David O'Donnell 2006
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <getopt.h>

unsigned NUMALLELES = 2;
bool be_quiet = false;

using namespace::std;

void PrintHelpText()
{
  cout <<
      "----------------------------------------------------------------------------"
      << endl <<
      "This program converts phased HapMap data to ADMIXMAP format" << endl
      << "Usage: convertdata -c <chrm> [-p <pop>] [-g <genotypesfile>]" << endl
      << " [-l <locusfile>] [-n <numloci>]" << endl
      << "where <pop> is the directory containing data files ('Eur', 'Afr' or 'Asian')," << endl
      << " defaults to current directory; " << endl
      << "chrm is a chromosome number. -c0 converts all chromosomes." << endl
      << "numloci is the number of loci per chromosome to use, default is all" << endl
      << "genotypesfile and locusfile default to 'genotypes.txt' and 'loci.txt'" << endl
      << "Note that HapMap input files must be named 'chrN_phased.txt', 'chrN_sample.txt'" << endl
      << "and 'chrN_legend.txt', where N is an integer from 1 to 22"
      << endl << endl;
}

int main(int argc, char **argv)
{
  if (argc < 3) {
    PrintHelpText();
    exit(0);
  }
  int ich;
  unsigned CHRNUM = 0;          //chromosome number
  string popprefix = ".";
  ofstream genotypesfile;
  ofstream locusfile;
  unsigned long userloci = 1000000000;  //some number > number of HapMap loci

  while ((ich = getopt(argc, argv, "hc:g:l:p:n:q")) != EOF) {
    switch (ich) {
    case 'h':{
        PrintHelpText();
      }
    case 'c':{
        CHRNUM = atoi(optarg);
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
    case 'n':{                 //number of loci
        const unsigned temp = atoi(optarg);
        if (temp > 0)
          userloci = temp;
        break;
      }
    case 'q':{                 //quiet mode
        be_quiet = true;
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

  //ofstream outcomefile("outcome_halfmissing.txt");
  string SNPID;
  double position = 0.0, prev = 0.0;
  string scrap;
  unsigned int locus = 0;
  // unsigned chromosome = 0;
  unsigned lastchr = 22;
  if (CHRNUM == 0)
    CHRNUM = 1;
  else
    lastchr = CHRNUM;
  vector < unsigned long >NUMLOCI;
  unsigned long TOTALLOCI = 0;

// ************************** PHASE 1: write locus file and header of genotypesfile, count number of loci  *********************
  genotypesfile << "\"Gameteid\"\t";
  locusfile << "\"SNPid\"\t\"NumAlleles\"\t\"DistanceinMb\"\n"; //locusfile header
  locusfile << setiosflags(ios::fixed) << setprecision(8);

  long unsigned abslocus = 0;   //absolute locus number
  ifstream legendfile;
  for (unsigned chr = CHRNUM; chr <= lastchr; ++chr) {
    cout << "\nChromosome " << chr << "  " << endl;
    locus = 0;
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

    getline(legendfile, scrap); //skip header
    while ((!legendfile.eof())) {
      prev = position;
      cout << "\rLocus    " << locus + 1 << flush;
      //we want cols: 0(snpid), 1(position in basepairs)
      legendfile >> SNPID >> position;
      // unsigned indiv_index = 0;
      if (SNPID.find_first_not_of(" \t\n\r") != string::npos) { //check for empty lines

        //write locusfile
        if ((position - prev) > 0) {    //strictly greater than to avoid having comp loci
          if (locus < userloci) {
            locusfile << SNPID << "\t" << 2 << "\t";
            if (locus == 0)
              locusfile << "#"; //missing value for first distance on chr
            else
//converts position in basepairs to distance in centiMorgans 
              //locusfile << 0.0000013 * (position-prev);
//convert position to distance in megabases
              locusfile << 0.000001 * (position - prev);
            locusfile << /*chrnumber << */ endl;
            //write header of genotypesfile
            genotypesfile << SNPID << "\t";
          }
          ++locus;
        }
//        else{//locus out of sequence
//          ++NUMBADLOCIONCHR;
//          ++NUMBADLOCI;
//          BadLoci.push_back(true);//excludes loci out of sequence
//          badlocusfile << SNPID << "\t" << chr << endl;
//      }
        ++abslocus;
      }
      SNPID.clear();
      //skip rest of line
      getline(legendfile, scrap);
    }
    legendfile.close();
    cout << endl << locus << " Loci ";
    NUMLOCI.push_back(locus);   //count number of loci
    TOTALLOCI += locus;
  }                             //end chr loop
  locusfile.close();
  cout << "\nFinished writing locusfile. " << TOTALLOCI << " loci" << endl;
  cout << endl;
  //}

  genotypesfile << endl;

// ************************** PHASE 2: Write genotypesfile  *********************


//Note: sample file is arranged with founders first, then children
//children are omitted from phased file
  string ID;
//       if(indiv%2) outcomefile << "#\n";
//       else outcomefile << 9 << endl;

  abslocus = 0;

//open first phased file to use to loop through gametes
//NOTE: assuming the same set of individuals for all chromosomes
  stringstream ss;
  ss << popprefix << "/chr" << CHRNUM << "_phased.txt";
  ifstream phasedfile(ss.str().c_str());
  if (!phasedfile.is_open()) {
    cout << "Could not open phased file " << ss.str() << endl;
    exit(1);
  }
  ss.clear();
  ss.str("");
  ss << popprefix << "/chr" << CHRNUM << "_sample.txt";
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
    if (!be_quiet)
      cout << "\nChromosome " << CHRNUM << "  " << flush;
    if (!(gamete % 2))
      samplefile >> ID >> scrap;
    else
      ID.append("_2");
    if (!be_quiet)
      cout << "\n" << gamete + 1 << " " << ID << " " << flush;
    //write indiv id and sex
    genotypesfile << ID << "\t";        // << sex[indiv] <<"\t";
//write genotypes for first chromosome
    for (unsigned j = 0; j < NUMLOCI[0]; ++j) {
      if (j < userloci)
        genotypesfile << allele + 1 << " ";
      phasedfile >> allele;
    }
//now loop through other chromosomes 
    for (unsigned chr = CHRNUM + 1; chr <= lastchr; ++chr) {
      if (!be_quiet)
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
      for (unsigned j = 0; j < NUMLOCI[chr - CHRNUM]; ++j) {
        phasedfile2 >> allele;
        if (j < userloci)
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

////////////////////////////////////////////////////////////////////////////
//  Program to convert ADMIXMAP input files to ANCESTRYMAP format.
//  Copyright(c) David O'Donnell 2005
//  
// 
// This code may be used in compiled form in any way you desire. This
// file may be redistributed and/or modified PROVIDING it is
// not sold for profit without the author's written consent, and
// providing that this notice and the author's name is included.
//
// This file is provided 'as is' with no expressed or implied warranty.
// The author accepts no liability if it causes any damage to your computer.
//
////////////////////////////////////////////////////////////////////////////

/*
  IMPORTANT: loci must all be SNPs, there must be only a single binary outcome variable
  and 2 populations.
  These requirements are NOT checked
  All input files must end with an end-of-line character for program to work properly.
  Set GETDIMS to get round this.
*/

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <stdexcept>

//#define GETDIMS 1
//set this to prompt for num individuals and num loci (assumed the same for all files)
//unset to determine from files (numloci allowed to be different)

#define CHECKIND 1
//set this to check num individuals in outcome file and genotypesfile match
//if they do not, program exits with error.
//  unset to determine num individuals only from outcome file

#define CHECKLOCI 1
//set this to check num loci in genotypesfile and locusfile match
// if not program issues warning but continues.
// unset to use value determined from genotypesfile, if specified.
//NB: program assumes numloci and locusnames in locusfile and historicallelefreqfile match
//    should check this too


#define CHECKLOCINAMES 1
//set this to check if locus names in genotypesfile and locus file match
//if not, program exits with error.

using namespace::std;

string dequote(const string& str)
//removes quotes from strings
{
  if (str.size() >= 2 && (str[0] == '"' || str[0]== '\'') && (str[str.size() - 1] == '"' || (str[str.size() - 1] == '\''))) {
    return str.substr(1, str.size() - 2);
  }

  return str;
}

string genotype2allelecount(const string& str){
  //converts genotype strings to allele counts (as strings)
  string astr;
  if(str == "1,1")astr = "0";
  else if(str == "2,2") astr = "2";
  else if(str == "") astr = "-1" ;
  else astr = "1";
  return astr;
}

static void openFile(ifstream *in, const char *fname)
{
  if (0 == fname || 0 == strlen(fname)) return;
  try{
    in->open(fname);
    if (!in->is_open()) {
      string msg = "Cannot open file for reading: \"";
      msg += fname;
      msg += "\"";
      throw (msg);
    }
  }
  catch(string msg){
    cout << msg <<endl;
    exit(1);
  }
}

static void openFile(ofstream *in, const char *fname)
{
  if (0 == fname || 0 == strlen(fname)) return;
  try{
    in->open(fname);
    if (!in->is_open()) {
      string msg = "Cannot open file for writing: \"";
      msg += fname;
      msg += "\"";
      throw (msg);
    }
  }
  catch(string msg){
    cout << msg <<endl;
    exit(1);
  }
}


int main()
{
  cout << "\t***************************************************" <<endl;
  cout << "\t ADMIXMAP to ANCESTRYMAP format file converter" << endl;
  cout << "\tCopyright (c) David O'Donnell 2005 " << endl;
  cout << "\t***************************************************" <<endl << endl;

  cout << "IMPORTANT: loci must all be SNPs, there must be only a single binary outcome variable"
       << " and 2 populations." <<endl
       << "All input files must end with an end-of-line character." << endl;

  //should determine these from data files
  //numInd is #rows of outcomevarfile or genotypesfile
  //numLoci is #cols of genotypesfile or #rows locusfile or 0.5*#rows allelefreqfile
  int numInd, numInd2, numLoci, numLoci2;
  string filename;
  string inputdir;
  string targetdir;
  bool dogeno = false, doind = false;
  string scrap; //for discarded info
  string *LociNames;

  cout << "Enter pathname of ADMIXMAP data directory"<<endl;
  cin >>  inputdir;
  inputdir = "./" + inputdir;
  cout << "Enter pathname of ANCESTRYMAP data directory"<<endl;
  cin >>  targetdir;
  targetdir = "./" + targetdir;

  /***************************************************
    Convert genotypesfile and outcomevarfile
  ***************************************************/

  //input files
  ifstream genotypesfile, outcomefile;
  //output files
  ofstream anc_genotypesfile, indivfile;

  cout << "Enter name of ADMIXMAP genotypes file or 'n' to skip"<<endl;
  cin >> filename;
  if(filename != "n") dogeno = true;
  if(dogeno){
    openFile(&genotypesfile, (inputdir + "/" + filename).c_str());
    cout << "Enter name of ANCESTRYMAP genotypes file"<<endl;
    cin >> filename;
    openFile(&anc_genotypesfile, (targetdir + "/" + filename).c_str());
  }

  cout << "Enter name of ADMIXMAP outcomevar file or 'n' to skip"<<endl;
  cin >> filename;
  if(filename != "n") doind = true;
  if(doind){
    openFile(&outcomefile, (inputdir + "/" + filename).c_str());
    cout << "Enter name of ANCESTRYMAP individual file"<<endl;
    cin >> filename;
    openFile(&indivfile, (targetdir + "/" + filename).c_str());
  }
  cout << endl;

#ifndef GETDIMS
  if(doind){
    //determine number of individuals by counting lines in outcomevar file
    numInd = -2; //to account for header and final newline
    while (! outcomefile.eof() )
      {
	getline (outcomefile, scrap);
	numInd++;
      }
    //rewind to start of file
    outcomefile.clear(); // to clear fail status caused by eof
    outcomefile.seekg(0);
    cout << numInd << " Individuals in outcome file" <<endl;
  }

  if(dogeno){
    //determine number of loci by counting quote characters in genotypes header
    numLoci = 0;
    char c, endline = '\n', quote = '\"';
    do{
      c = genotypesfile.get();
      if(c == quote)++numLoci;
    }
    while(c != endline);
    numLoci = numLoci/2 -1;
    if(numLoci < 1) {
      cout << "header strings must be enclosed in double quotes" << endl;
      exit(1);
    }
    cout << numLoci << " loci in genotypes file" <<endl;

    //determine number of individuals if no outcome file
#ifndef CHECKIND
    if(!doind){
#endif
      numInd2 = -1; //no header this time, already read
      while (! genotypesfile.eof() )
	{
	  getline (genotypesfile, scrap);
	  numInd2++;
	}
    cout << numInd2 << " Individuals in genotypes file" <<endl;
    if(doind){
      if(numInd != numInd2){
	cerr << "ERROR: number of individuals in outcomevar file and genotypes file do not match"<<endl;
	exit(1);
      }
    }
    else numInd = numInd2;
#ifndef CHECKIND
    }
#endif
    //rewind to start of file
    genotypesfile.clear();
    genotypesfile.seekg(0);
  }
#endif

  if(doind || dogeno){
#ifdef GETDIMS
    cout << "Enter number of Individuals: ";
    cin >> numInd;
    if(dogeno){
      cout << "Enter number of Loci: ";
      cin >> numLoci;
    }
#endif
    cout << "Processing ..." << endl;
    
    int outcomevar;
    string Indiv_ID;
    string **genotypes;
    string allele_string;
    
    //write header lines
    if(doind){
      indivfile << "Indiv_ID\tGender\tStatus" << endl;
    }
    if(dogeno){
      anc_genotypesfile << "SNP_ID\tIndiv_ID\tVart_allele_cnt" << endl;
    }
    
    if(doind){
      //discard header of outcomevar file
      getline(outcomefile,scrap);
    }
    
    if(dogeno){
      //read header of genotypesfile
      LociNames = new string[numLoci+1];
      //NB LociNames[0] = "IDs"
      for(int i = 0; i <= numLoci; ++i){
	genotypesfile >> LociNames[i];
      }
      //allocate genotypes array
      genotypes = new string *[numInd*numLoci];
      for(int s = 0; s < numInd*numLoci; ++s) genotypes[s] = new string[3];
      
    }
    
    
    for(int i = 0; i < numInd; ++i){
      if(dogeno){
	//read indiv ID
	genotypesfile >> Indiv_ID;
      }
      if(doind){
	if(!dogeno){
	  //write Indiv and U(nknown) Gender to Indivfile
	  indivfile << "Ind_"<<i+1 <<"\t"<<"U\t";
	}
	else{
	  //write Indiv and U(nknown) Gender to Indivfile
	  indivfile << dequote(Indiv_ID) <<"\t"<<"U\t";
	}
	//write Status to Indivfile
	outcomefile >> outcomevar;
	indivfile << (outcomevar==0?"Control":"Case") <<endl;
      }
      if(dogeno){
	//process each line of genotypesfile
	for(int j = 0; j < numLoci; ++j){
	  genotypesfile >> allele_string;
	  genotypes[numInd*j +i][0] = dequote(LociNames[j+1]);//locus name
	  genotypes[numInd*j +i][1] = dequote(Indiv_ID);      //Indiv name
	  genotypes[numInd*j +i][2] = genotype2allelecount(dequote(allele_string));//variantallelecount
	  //NB: using allele2 as variant allele here
	}
      }
    }
    
    if(dogeno){
      //write genotype data
      for(int i = 0; i < numInd*numLoci; ++i)
	anc_genotypesfile << genotypes[i][0] << "\t" << genotypes[i][1] << "\t" << genotypes[i][2] <<endl;
      
      //clean up
      genotypesfile.close();
      anc_genotypesfile.close();
      for(int s = 0; s < numInd*numLoci; ++s) delete[] genotypes[s];
      delete[] genotypes;
    }
    
    if(doind){
      outcomefile.close();
      indivfile.close();
    }
  }
  
  /********************************************************************************
    Convert locusfile and allelefreqfile
  ********************************************************************************/
  bool domarkers = false;
  bool isprior = false;
  //input files
  ifstream locusfile, allelefreqfile;
  //outputfile
  ofstream markerfile;


  cout << "Enter name of locus file or 'n' to skip"<<endl;
  cin >> filename;
  domarkers = (filename != "n");
  string filename2;
  if(domarkers){
    openFile(&locusfile, (inputdir + "/" + filename).c_str());
    cout << "Enter name of historicallelefreqs or 'n' to skip"<<endl;
    cin >> filename2;
    if(filename2 == "n"){
      cout << "Enter name of priorallelefreqs or 'n' to skip"<<endl;
      cout<< "Note that priors are assumed to be counts+0.5\n";
      cin >> filename2;
      if(filename2 == "n") domarkers = false;
      else isprior = true;
      }
  }

  if(domarkers){
    openFile(&allelefreqfile, (inputdir + "/" + filename2).c_str());
    cout << "Enter name of ANCESTRYMAP marker file"<<endl;
    cin >> filename;
    openFile(&markerfile, (targetdir + "/" + filename).c_str());

#ifndef GETDIMS
#ifndef CHECKLOCI
    if(!dogeno){
#endif
      //determine number of loci (if not already known from genotypesfile) by counting lines of locusfile
      //should check allelefreqfile too
      numLoci2 = -2; //to account for header and final newline
      while (! locusfile.eof() )
	{
	  getline (locusfile, scrap);
	  numLoci2++;
	}
    //rewind to start of file
      locusfile.clear();
      locusfile.seekg(0);
      cout << numLoci2 << " Loci in locus file" <<endl;
      if(dogeno){
	if(numLoci2 != numLoci){
	  cerr << "WARNING: Numbers of loci in genotypes file and locusfile do not match"<<endl;
	  //exit(1);
	}
      }
#ifndef CHECKLOCI
      }
    else numLoci2 = numLoci;
#endif
#endif

#ifdef GETDIMS
    if(!dogeno){
      cout << "Enter Number of Loci: ";
      cin >> numLoci2;
    }
    else numLoci2 = numLoci;
#endif

    //should really check allelefreqfile too

    cout << endl << "Creating markerfile..." << endl;

    //discard headers
    getline(locusfile, scrap);
    getline(allelefreqfile, scrap);

    double mapdist;    //map distance
    double GenPos;    //genetic position
    int PhysPos = 0; //dummy for physical position
    int chr = 0;      //chromosome number
    string snpid, snpid2;     // SNP label
    double PopA1, PopA2, PopB1, PopB2; //allele counts

    //write header line
    markerfile << "SNP_ID \t Chr \t Gen_Pos \t Phys_Pos \t  PopA_vart_cnt \t PopA_ref_cnt \t PopB_vart_cnt \t PopB_ref_cnt"<<endl;

    for(int i = 0; i < numLoci2; ++i){
      //process line of locusfile
      locusfile >> snpid;

#ifdef CHECKLOCINAMES
      if(dogeno){
	//check locus names
	if(snpid != LociNames[i+1]) {
	  cerr << "Locus Names in genotypes file and locusfile do not match!" <<endl;
	  exit(1);
	}
      }
#endif

      //discard numalleles (always 2)
      locusfile >> scrap;
      //read map distance and determine genetic position (on chromosome)
      locusfile >> mapdist;
      if(mapdist == 100){ //not ideal but will have to do. Perhaps read as string, then test, then convert to double.
	++chr;//new chromosome
	GenPos = 0.0;
      }
      else GenPos += mapdist;
      markerfile <<  dequote(snpid) <<"\t" << chr <<"\t";
      markerfile.precision(10);
      markerfile << GenPos;
      markerfile.precision(6);
      markerfile <<"\t" << ++PhysPos <<"\t";

      //read 2 lines of allelefreqfile
      //NB: using allele1 as "reference" and allele2 as "variant"
      allelefreqfile >> snpid >> PopA2 >> PopB2 >> snpid2 >> PopA1 >> PopB1;
      if(isprior) {PopA2 -= 0.5; PopB2 -= 0.5; PopA1 -= 0.5; PopB1 -= 0.5;}

      //check locus names
      //if(snpid != snpid2) {cout << "Error in locusfile: loci names should be in pairs"<<endl;}//break;}
      //if(snpid != LociNames[i+1])
      //  {cout<<"Loci Names in genotypes file and allelefreqfile do not match!"<<endl;}//break;}
      //should check match with locus file too

      //write allelecounts
      markerfile << (int)PopA1 <<"\t" << (int)PopA2 <<"\t"<< (int)PopB1 <<"\t"<< (int)PopB2 <<endl;
    }

    //clean up
    locusfile.close();
    allelefreqfile.close();
    markerfile.close();
  }

  cout << "\nConversion complete" <<endl;
  cout << "Please check headers of output files"<<endl;

  //end program
  return 0;
}

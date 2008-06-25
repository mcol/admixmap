////////////////////////////////////////////////////////////////////////////
//  Program to convert ANCESTRYMAP input files to ADMIXMAP format.
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
  IMPORTANT:
  All input files must end with an end-of-line character for program to work properly.
  Set GETDIMS to get round this.
  Genotypes file is assumed to be sorted by locus, then by individual as per manual.
  All ANCESTRYMAP strings are assumed to be WITHOUT QUOTES. If this is not the case, use the dequote function.
  There is no check that locus labels in genotypes file and marker file match or that individual IDs in genotypes
  file and individual files match (but these are easy to do).
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
//only meaningful if GETDIMS is undef'd
//set this to check num individuals in individual file and genotypesfile match
//if they do not, program exits with error.
//  unset to determine num individuals only from outcome file

#define CHECKLOCI 1
//only meaningful if GETDIMS is undef'd
//set this to check num loci in genotypesfile and markerfile match
// if not program issues warning but continues.
// unset to use value determined from genotypesfile, if specified.

using namespace::std;

//removes quotes from strings
string dequote(const string& str)
{
  if (str.size() >= 2 && (str[0] == '"' || str[0]== '\'') && (str[str.size() - 1] == '"' || (str[str.size() - 1] == '\''))) {
    return str.substr(1, str.size() - 2);
  }

  return str;
}

//converts allelecounts in ANCESTRYMAP genotypes file to genotype strings
string allelecount2genotype(const string& str){
  string astr;
  if(str == "0")astr = "\"1,1\"";
  else if(str == "2") astr = "\"2,2\"";
  else if(str == "1") astr = "\"1,2\"" ;
  else astr = "\"\"";//missing genotype
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
  cout << "\tANCESTRYMAP- to ADMIXMAP-format file converter" << endl;
  cout << "\tCopyright (c) David O'Donnell 2005 " << endl;
  cout << "\t***************************************************" <<endl << endl;
  cout << "IMPORTANT: All input files must end with an end-of-line character for program to work properly." << endl;

  int numInd, numInd2, numLoci;
  string filename;
  string inputdir;
  string targetdir;
  string scrap;//for discarded info
  bool dogeno = false, doind = false;

  cout << "Enter pathname of ANCESTRYMAP data directory"<<endl;
  cin >>  inputdir;
  inputdir = "./" + inputdir;
  cout << "Enter pathname of ADMIXMAP data directory"<<endl;
  cin >>  targetdir;
  targetdir = "./" + targetdir;

  //input files
  ifstream anc_genotypesfile, indivfile;
  //output files
  ofstream genotypesfile, outcomefile;


  cout << "Enter name of ANCESTRYMAP genotypes file or 'n' to skip"<<endl;
  cin >> filename;
  if(filename != "n") dogeno = true;
  if(dogeno){
    openFile(&anc_genotypesfile, (inputdir + "/" + filename).c_str());
    cout << "Enter name of ADMIXMAP genotypes file"<<endl;
    cin >> filename;
    openFile(&genotypesfile, (targetdir + "/" + filename).c_str());
  }

  cout << "Enter name of ANCESTRYMAP individual file or 'n' to skip"<<endl;
  cin >> filename;
  if(filename != "n") doind = true;
  if(doind){
    openFile(&indivfile, (inputdir + "/" + filename).c_str());
    cout << "Enter name of ADMIXMAP outcomevar file"<<endl;
    cin >> filename;
    openFile(&outcomefile, (targetdir + "/" + filename).c_str());
  }
  cout << endl;

#ifndef GETDIMS
  if(doind){
    //determine number of individuals by counting lines in individual file
    numInd = -2; //to account for header and final newline
    while (! indivfile.eof() )
      {
	getline (indivfile, scrap);
	numInd++;
      }
    //rewind to start of file
    indivfile.clear(); // to clear fail status caused by eof
    indivfile.seekg(0);
    cout << numInd << " Individuals in individual file" <<endl;
  }

  if(dogeno){
    //determine number of loci and number of individuals in genotypesfile
    getline(anc_genotypesfile,scrap);//header
    streampos line1 = anc_genotypesfile.tellg();//save position
    string label1, label;
    int numlines = 0;
    numLoci = 0;
    anc_genotypesfile >> label1;
    //determine num individuals by counting lines with first locus repeated
    do{
	getline(anc_genotypesfile, scrap);
	anc_genotypesfile >> label;
	++numlines;
      }
    while(label == label1);
    numInd2 = numlines;
    while(! anc_genotypesfile.eof())
      {
	getline(anc_genotypesfile, scrap);
	++numlines;
      }
    numLoci = (numlines - 1)/numInd2;
    cout << numLoci<< " loci in genotypesfile"<<endl;
 
    //determine number of individuals if no outcome file
#ifndef CHECKIND
    if(!doind){
#endif
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
    anc_genotypesfile.clear();
    anc_genotypesfile.seekg(line1);
  }
#endif

  if(dogeno || doind){
#ifdef GETDIMS
    cout << "Enter number of Individuals: ";
    cin >> numInd;
    if(dogeno){
      cout << "Enter number of Loci: ";
      cin >> numLoci;
      //discard header of genotypes file
      getline(anc_genotypesfile,scrap);
    }
#endif
    if(doind){
      cout << "Enter name of outcome variable:" <<endl;
      cin >> scrap ;
      //outcomefile << "\"Outcome\"" << endl;
      //write header of outcome file
      outcomefile << "\"" << scrap << "\"" << endl;
      //discard header of outcomevar file
      getline(indivfile,scrap);
    }
    cout <<endl<< "Processing ..." << endl;
  }

  if(dogeno){
 
    //write "ID" to genotypes file
    genotypesfile << "\"IDs\" ";
  }

  string outcomevar;
  string Indiv_ID;
  string ***genotypes;
  string allele_string;


  /*
    Convert genotypesfile and individual file
  */
  if(dogeno){
    //allocate genotypes array
    genotypes = new string**[numInd];
    for(int i = 0; i < numInd;++i){
      genotypes[i] = new string*[numLoci];
      for(int j = 0; j < numLoci; ++j)genotypes[i][j] = new string[3];
    }
    //read genotype data
    for(int j = 0; j < numLoci; ++j){
      for(int i = 0; i < numInd; ++i){
	anc_genotypesfile >> genotypes[i][j][0] >> genotypes[i][j][1];

	//determine if missing genotype by checking for non-whitespace chars
	getline(anc_genotypesfile, scrap);
	int pos = scrap.find_first_not_of(" \t\r\n");
	if(pos==string::npos){//missing
	  allele_string = "";
	}
	else{//not missing
	  allele_string = scrap.substr(pos, scrap.find_first_of(" \t\r\n", pos));
	}

	genotypes[i][j][2] = allelecount2genotype(allele_string);
      }
	//write locus names to genotypes file header
	genotypesfile << "\"" << genotypes[0][j][0] << "\" ";
    }
    genotypesfile <<endl;
  }

  for(int i = 0; i < numInd; ++i){
    if(dogeno){
      //write indiv ID
      genotypesfile << "\"" << genotypes[i][0][1] << "\"\t";
      //write each line of genotypesfile
      for(int j = 0; j < numLoci; ++j){
	genotypesfile << genotypes[i][j][2] << " ";
      }
    
      genotypesfile <<endl;
    }

    if(doind){
      //read line of individual file
      indivfile >> Indiv_ID >> scrap >> outcomevar;
      //write outcome to outcomevar file
      if(outcomevar == "Control") outcomefile << 0 << endl;
      else if(outcomevar == "Case") outcomefile << 1 << endl;
      else {cout << "Unexpected status found in individual file: " << outcomevar <<endl;
	outcomefile << outcomevar <<endl;}
    }
  }

  //clean up
  if(dogeno){
    genotypesfile.close();
    anc_genotypesfile.close();
    for(int i = 0; i < numInd;++i){
      for(int j = 0; j < numLoci; ++j)delete[]genotypes[i][j];
      delete[] genotypes[i];
    }
    delete[] genotypes;
  }

  if(doind){
    outcomefile.close();
    indivfile.close();
  }

  /*
    Convert markerfile
  */
  int numLoci2;
  bool domarkers = false;
  bool dopriors = false; bool dohistoric = false;
  //input files
  ofstream locusfile, allelefreqfile, priorallelefreqfile;
  //outputfile
  ifstream markerfile;

  cout << "Enter name of ANCESTRYMAP marker file or 'n' to skip"<<endl;
  cin >> filename;
  domarkers = (filename != "n");

  if(domarkers){
    openFile(&markerfile, (inputdir + "/" + filename).c_str());

#ifndef GETDIMS
#ifndef CHECKLOCI
    if(!dogeno){
#endif
      //determine number of loci (if not already known from genotypesfile) by counting lines of markerfile
      numLoci2 = -2; //to account for header and final newline
      while (! markerfile.eof() )
	{
	  getline (markerfile, scrap);
	  numLoci2++;
	}
    //rewind to start of file
      markerfile.clear();//to clear fail status caused by eof
      markerfile.seekg(0);
      cout << numLoci2 << " Loci in marker file" <<endl;
      if(dogeno){
	if(numLoci2 != numLoci){
	  cerr << "WARNING: Numbers of loci in genotypes file and markerfile do not match"<<endl;
	  //exit(1);
	}
      }
#ifndef CHECKLOCI
      }
    else numLoci2 = numLoci;
#endif
#endif

    cout << "Enter name of ADMIXMAP locus file"<<endl;
    cin >> filename;
    openFile(&locusfile, (targetdir + "/" + filename).c_str());
    cout << "Enter name of ADMIXMAP historicallelefreqs file or 'n' to skip"<<endl;
    cin >> filename;
    dohistoric = (filename != "n");
    openFile(&allelefreqfile, (targetdir + "/" + filename).c_str());
    cout << "Enter name of ADMIXMAP priorallelefreqs file or 'n' to skip"<<endl;
    cin >> filename;
    dopriors = (filename != "n");
    openFile(&priorallelefreqfile, (targetdir + "/" + filename).c_str());


#ifdef GETDIMS
    if(!dogeno){
      cout << "Enter Number of Loci: ";
      cin >> numLoci2;
    }
    else numLoci2 = numLoci;
#endif

    cout << endl << "Converting markerfile..." << endl;
    //discard headers
    getline(markerfile, scrap);

    double mapdist;    //map distance
    double GenPos, prevGenPos;    //genetic position
    int PhysPos = 0; //dummy for physical position
    int chr = 0, prevchr = 0;      //chromosome number
    string snpid;     // SNP label
    int PopA1, PopA2, PopB1, PopB2; //allele counts

    //prompt for popnames (optional)
//     string p1, p2;
//     cout << "Enter name of Population 1: " ;
//     cin >> p1;
//     cout << "Enter name of Population 2: " ;
//     cin >> p2;

    //write header lines
    locusfile << "\"SNP_ID\" \t \"NumAlleles\" \t \"MapDist\"" << endl;
    if(dohistoric)allelefreqfile << "\"SNP_ID\" \t \"PopA\" \t \"PopB\"" << endl;
    if(dopriors)priorallelefreqfile << "\"SNP_ID\" \t \"PopA\" \t \"PopB\"" << endl;
    //allelefreqfile << "\"SNP_ID\" \t \""<<p1<<"\" \t \""<<p2<<"\"" << endl;

    for(int i = 0; i < numLoci2; ++i){
      //process line of markerfile
      markerfile >> snpid >> chr >> GenPos >> PhysPos >> PopA1 >> PopA2 >> PopB1 >> PopB2;

      //check locus names
      //if(snpid != LociNames[i+1])
      //  {cout<<"Loci Names in genotypes file and markerfile do not match!"<<endl;}//break;}
      locusfile << "\"" << snpid << "\"" << "\t" << 2 << "\t";

      //read genetic position and determine map distance
      if(chr == prevchr){
	mapdist = GenPos - prevGenPos;
      }
      else{
	mapdist = 100; //new chromosome
      }

      locusfile << mapdist << endl;

      //write 2 lines of allelefreqfile
      //NB: using allele1 as "reference" and allele2 as "variant"
      if(dohistoric){
        allelefreqfile << "\"" << snpid << "\"" << "\t" << PopA2 << "\t" << PopB2 << endl;
        allelefreqfile << "\"" << snpid << "\"" << "\t" << PopA1 << "\t" << PopB1 << endl;
        }
      if(dopriors){
        priorallelefreqfile << "\"" << snpid << "\"" << "\t" << PopA2+0.5 << "\t" << PopB2+0.5 << endl;
        priorallelefreqfile << "\"" << snpid << "\"" << "\t" << PopA1+0.5 << "\t" << PopB1+0.5 << endl;
        }

      prevGenPos = GenPos;
      prevchr = chr;
    }

    //clean up
    locusfile.close();
    allelefreqfile.close();
    priorallelefreqfile.close();
    markerfile.close();
  }

  cout << "Conversion complete" <<endl;
  cout << "Please check headers of output files"<<endl;

  //end program
  return 0;
}

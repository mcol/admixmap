/*
prog to convert hapmap data (for a single chromsome) to ADMIXMAP format
supply chromosome number as arg. 0 (default) means all.
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

unsigned CHRNUM = 0;//chromosome number
unsigned NUMIND = 90;
unsigned NUMVALIDINDS =0;
unsigned long NUMLOCI = 0;
unsigned long NUMBADLOCI = 0;
unsigned long NUMBADLOCIONCHR = 0;
unsigned NUMALLELES = 2;

using namespace::std;
string getGenotype(string g, string a){
  //converts a base pair string g, eg "CT", to an admixmap format genotype eg "1,2"
  //a is the reference string of alleles eg "C/T"
  stringstream s;
  s << "\"" << (int)(g[0]==a[0])+2*(int)(g[0]==a[2]) << "," << (int)(g[1]==a[0])+2*(int)(g[1]==a[2]) << "\"";
  return s.str();
}

int main(int argc, char **argv){
    if(argc==1){
	cout << "Usage: convertdata c [genotypesfile [locusfile]]" << endl
	     << "where c is a chromosome number; chr=0 converts all chromosomes." << endl;
	exit(0);
    }
  if(argc>1) CHRNUM = atoi(argv[1]);
  ifstream genotypesin;
  ifstream infofile("pedinfo2sample_CEU.txt");//custom file with info on individuals' sex and validity (only unrelated inds are valid)
  ofstream genotypesfile;
  if(argc>2)genotypesfile.open(argv[2], ios::out);
  else genotypesfile.open("genotypes.txt", ios::out);
  ofstream locusfile;
  if(argc>3)locusfile.open(argv[3]);
  else locusfile.open("loci.txt");

  ofstream badlocusfile("excludedSNPs.txt");
  //ofstream outcomefile("outcome_halfmissing.txt");
  string alleles;
  string genotype;
  string obs;//observed genotype
  string SNPID;
  vector<string> BadLoci;
  vector<string> INDIVID;
  vector<int> indisvalid;
  vector<int> sex;//1 = male, 2 = female
  double position = 0.0, prev = 0.0;
  string scrap;
  int locus = 0;
  int indiv = 0;
  unsigned chromosome = 0;
  unsigned lastchr = 22;
  if(CHRNUM==0)CHRNUM = 1;
  else lastchr = CHRNUM;

  //NB: assumes individuals in same order in all files

  //read info file for indiv id, sex
  NUMIND = 0;

  while(!infofile.eof()){
    int _sex, valid;
    string _id;
    infofile >> _id >> _sex >> valid;
    if(_id.find_first_not_of(" \t\n\r")!=string::npos){//check for empty lines
      INDIVID.push_back(_id);
      sex.push_back(_sex);
      indisvalid.push_back(valid);
      //cout << _id << " " << valid << endl;
      ++indiv;
      if(valid)++NUMVALIDINDS;
    }
  }
  NUMIND = indiv;
  infofile.close();
  cout << NUMIND << " individuals" <<endl;
  cout << NUMVALIDINDS << " unrelated individuals" << endl;

  // *** write locus file and header of genotypesfile, count number of loci 
  genotypesfile << "\"Individ\"\t";
  locusfile << "\"SNPid\"\t\"NumAlleles\"\t\"Distance\"\n";
  locusfile << setiosflags(ios::fixed) << setprecision(8);

  for(unsigned chr = CHRNUM; chr <= lastchr; ++chr){
      cout << "\nChromosome " << chr << "  " <<endl;
      NUMBADLOCIONCHR = 0;
      locus = 0;
    stringstream ss;
    //ss << "gzip -d chr" << chr << ".gz";
    //cout << "\nUnzipping..."<<flush;
    //system(ss.str().c_str());
    //ss.clear();
    ss << "chr" << chr;
    genotypesin.clear();//to clear fail status at eof
    genotypesin.open(ss.str().c_str());
    //cout << "\nReading raw genotypes...\n"<<flush;
    getline(genotypesin, scrap);//skip header
    while (! genotypesin.eof() ){
      ////use this to split up large locusfiles
      //if(!(locus%10000)){
      //locusfile.close();
      //stringstream s; s << "SNP" << locus +10000  << ".txt";
      //locusfile.open(s.str().c_str());
      //}
      prev = position;
      cout << "\rLocus    " << locus+1 << flush;
      //we want cols: 0(snpid), 1(alleles), 3(position in basepairs), 11 onwards(indiv genotypes)
      genotypesin >> SNPID >> alleles
		  >> scrap >> position;
      if(SNPID.find_first_not_of(" \t\n\r")!=string::npos){//check for empty lines
	//write locusfile
	  if((position-prev)<0){
	      ++NUMBADLOCIONCHR;
	      ++NUMBADLOCI;//BadLoci.push_back(SNPID);//excludes loci out of sequence
	      badlocusfile << SNPID << "\t" << chr << endl;
	  }
	else{
	  locusfile << SNPID << "\t" <<  2 << "\t";
	  if(locus==0) locusfile << "#" ;
	  else 
	    locusfile << 0.0000013 * (position-prev);//converts position in basepairs to distance in Morgans 
	  locusfile << /*chrnumber <<*/ endl;
	  //write header of genotypesfile
	  genotypesfile << SNPID << "\t";
	  ++locus;
	}
	
      }
      SNPID.clear();       
      //skip rest of line
      getline(genotypesin, scrap);
    }
    genotypesin.close();
    cout << endl << locus << " Loci, " << NUMBADLOCIONCHR  << " loci excluded";
    NUMLOCI += locus;
  }
   locusfile.close();
   cout << "\nFinished writing locusfile. " << NUMLOCI <<  " loci" << endl;
   //if(BadLoci.size()){
       cout << NUMBADLOCI/*BadLoci.size()*/ << " loci excluded: ";
     //for(vector<string>::const_iterator i = BadLoci.begin(); i < BadLoci.end(); ++i)cout << *i << " "; 
     cout << endl;
     //}
     badlocusfile.close();

   genotypesfile << endl;
   // *** Write genotypesfile
   
    for(indiv = 0; indiv < NUMIND; ++indiv) if(indisvalid[indiv]){
      cout << "\n" << indiv+1 << " " << INDIVID[indiv] << " "  << flush;
      //write indiv id and sex
      genotypesfile << INDIVID[indiv] << "\t";// << sex[indiv] <<"\t";
      position = 0.0;     
//       if(indiv%2) outcomefile << "#\n";
//       else outcomefile << 9 << endl;
      for(unsigned chr = CHRNUM; chr <= lastchr; ++chr){
	  cout << "\nChromosome " << chr << "  "<< flush ;
	  stringstream ss;
	  ss << "chr" << chr;
	  genotypesin.open(ss.str().c_str());
 
	  //open genotypes input file
	  genotypesin.clear(); // to clear fail status caused by eof
	  //genotypesin.seekg(0);
	  getline(genotypesin, scrap);//skip header
      
	  //read genotypes, one locus (row) at a time
	  //for(locus = 0; locus < NUMLOCI+ NUMBADLOCI/*BadLoci.size()*/; ++locus){
	  prev = position;
	  genotypesin >> SNPID >> alleles >> scrap >> position;
	  while (! genotypesin.eof() ){
	      if(position-prev >=0.0){//skip loci out of sequence
		  //cout << "  locus" << locus+1 << " " << SNPID << " " <<position << " " << prev << endl;
		  for(int col = 0; col < 7+indiv; ++col)genotypesin >> scrap;//skip to col for this individual, genotypes start at col 11
		  genotypesin >> obs;
		  
		  //write to genotypesfile in admixmap format
		  //set 50% to missing at half the loci
// 		  if(indiv%2 && locus%2){
// 		      genotypesfile << "\"0,0\"" << "\t";
// 		  }
// 		  else{
 		      genotypesfile << getGenotype(obs, alleles) << "\t";
// 		  }
	      }
	      getline(genotypesin, scrap);//skip remaining individuals in line
	      prev = position;//start reading next line
	      genotypesin >> SNPID >> alleles >> scrap >> position;
	  }
	  genotypesin.close();
      }
      
      genotypesfile << endl;
    }
  genotypesfile.close();
//  outcomefile.close();
  cout << "\nFinished writing genotypesfile" << endl;
}

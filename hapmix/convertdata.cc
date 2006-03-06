/*
prog to convert hapmap data (for a single chromsome) to ADMIXMAP format
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

//#define CHRNUM 22;//chromosome number
unsigned NUMIND = 90;
unsigned NUMVALIDINDS =0;
unsigned NUMLOCI = 5;
unsigned NUMALLELES = 2;

using namespace::std;
string getGenotype(string g, string a){
  //converts a base pair string g, eg "CT", to an admixmap format genotype eg "1,2"
  //a is the reference string of alleles eg "C/T"
  stringstream s;
  s << "\"" << (int)(g[0]==a[0])+2*(int)(g[0]==a[2]) << "," << (int)(g[1]==a[0])+2*(int)(g[1]==a[2]) << "\"";
  return s.str();
}

int main(){
  ifstream genotypesin("genotypes_chr22_CEU.b35.txt");//edit chromosome number as appropriate
  ifstream infofile("pedinfo2sample_CEU.txt");//custom file with info on individuals' sex and validity (only unrelated inds are valid)
  ofstream genotypesfile("genotypes.txt");
  ofstream locusfile("loci.txt");

  string alleles;
  string genotype;
  string obs;//observed genotype
  string SNPID;
  vector<string> INDIVID;
  vector<int> indisvalid;
  vector<int> sex;//1 = male, 2 = female
  double position = 0.0, prev = 0.0;
  string scrap;
  int locus = 0;
  int indiv = 0;

  //NB: assumes individuals in same order in both files

  //read info file for indiv id, sex
  NUMIND = 0;
  while(!infofile.eof()){
    string _id;
    int _sex, valid;

    infofile >> _id >> _sex >> valid;
    INDIVID.push_back(_id);
    sex.push_back(_sex);
    indisvalid.push_back(valid);
    //cout << _id << " " << valid << endl;

    ++indiv;
  }
  NUMIND = indiv-1;
  infofile.close();
  cout << NUMIND << " individuals" <<endl;

  //write locus file and header of genotypesfile 
  genotypesfile << "\"Individ\"\t";
  locusfile << "\"SNPid\"\t\"NumAlleles\"\t\"Distance\"\n";
  locusfile << setiosflags(ios::fixed) << setprecision(8);

   getline(genotypesin, scrap);//skip header
   while (! genotypesin.eof() ){
     ////use this to split up large locusfiles
     //if(!(locus%10000)){
     //locusfile.close();
     //stringstream s; s << "SNP" << locus +10000  << ".txt";
     //locusfile.open(s.str().c_str());
     //}
     prev = position;
     cout << "\rLocus    " << locus << flush;
     //we want cols: 0(snpid), 1(alleles), 3(position in basepairs), 11 onwards(indiv genotypes)
     genotypesin >> SNPID >> alleles
 		>> scrap >> position;
     //write locusfile
     locusfile << SNPID << "\t" <<  2 << "\t";
     if(locus==0) locusfile << "#";
     else locusfile << 0.0000013 * (position-prev); 
     locusfile << /*chrnumber <<*/ endl;
     //write header of genotypesfile
     genotypesfile << SNPID << "\t";
 
     //skip rest of line
     getline(genotypesin, scrap);

     ++locus;
   }
   NUMLOCI  = locus-1;
   locusfile.close();
   cout << "\nFinished writing locusfile. " << NUMLOCI <<  " loci" << endl;

   genotypesfile << endl;
   for(indiv = 0; indiv < NUMIND; ++indiv) if(indisvalid[indiv]){
     ++NUMVALIDINDS;
     //cout << INDIVID[indiv] << " "  << indisvalid[indiv] << endl << flush;
     //write indiv id and sex
     genotypesfile << INDIVID[indiv] << "\t";// << sex[indiv] <<"\t";
     //rewind to start of file
     genotypesin.clear(); // to clear fail status caused by eof
     genotypesin.seekg(0);
     getline(genotypesin, scrap);//skip header
     
     //read genotypes, one locus (row) at a time
    for(locus = 0; locus < NUMLOCI; ++locus){
      genotypesin >> SNPID >> alleles;
      for(int col = 0; col < 9+indiv; ++col)genotypesin >> scrap;//skip to col for this individual, genotypes atart at col 11
      genotypesin >> obs;
      
      //write to genotypesfile in admixmap format
      genotypesfile << getGenotype(obs, alleles) << "\t";
      //cout << indiv << " "  << INDIVID[indiv] << " " << locus << " " << SNPID << " " << obs << " " <<  alleles << " " << getGenotype(obs, alleles) << endl << flush;
      //if(!(locus%10))system("pause");
      getline(genotypesin, scrap);//skip remaining individuals in line
    }
    genotypesfile << endl;
   }
   cout << "\nFinished writing genotypesfile" << endl;
   cout << NUMVALIDINDS << " unrelated individuals" << endl;
   
   genotypesfile.close();
   genotypesin.close();
   
}

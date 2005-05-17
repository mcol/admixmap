#include "Genome.h"
#include <stdlib.h>
#include "Chromosome.h"

using namespace std;

// composite-level methods

Genome::Genome()
{
  NumberOfCompositeLoci = 0;
  NumberOfChromosomes = 0;
  TheArray = 0;
  LengthOfGenome = 0;
  LengthOfXchrm = 0;
  TotalLoci = 0;
}

Genome::Genome( int size )
{
  if ( size < 1 ){
    size = 1;
  }
  NumberOfCompositeLoci = size;
  TheArray = new CompositeLocus*[ NumberOfCompositeLoci ];
  Distances.SetNumberOfElements( NumberOfCompositeLoci );
  LengthOfGenome = 0;
  LengthOfXchrm = 0;
  TotalLoci = 0;
}


Genome::~Genome()
{
  if(NumberOfCompositeLoci > 0){
    for(unsigned int i=0; i<NumberOfCompositeLoci; i++){
      //delete TheArray[i];
    }
    delete TheArray;
  }
  delete[] SizesOfChromosomes; 
}

//sets the labels of all composite loci in TheArray using SetLabel function in CompositeLocus
void Genome::SetLabels(const vector<string> &labels, Vector_d temp)
{
    int index = -1; // counts through number of composite loci
    int locus = 0;

    for (size_t count = 2; count < labels.size(); ++count) {
        const string& label = labels[count];
        if (temp.GetNumberOfElements() == 1 || temp(count - 1)) {
            index++;
            locus = 0;
            TheArray[index]->SetLabel(0,label);
        } else if( count != 0 ) {
            locus++;
            TheArray[index]->SetLabel(locus,label);
        }
    }
}

void Genome::SetLabel(int index, string str)
{
  for( unsigned int i = 0; i < NumberOfCompositeLoci; i++ ){
    TheArray[i]->SetLabel(index, str);
  }
}

//Creates an array of pointers to Chromosome objects
//also determines length of genome, NumberOfCompositeLoci, TotalLoci, NumberOfChromosomes, SizesOfChromosomes, 
//LengthOfXChrm and chrmandlocus
Chromosome** Genome::GetChromosomes( int populations, vector<string> chrmlabels)
{
  int *cstart = new int[NumberOfCompositeLoci];
  int *cfinish = new int[NumberOfCompositeLoci];
  int cnum = -1; //number of chromosomes -1
  int lnum=0;
  _chrmandlocus.resize(NumberOfCompositeLoci);
  X_data = false;
  
  TotalLoci = 0;
  for(unsigned int i=0;i<NumberOfCompositeLoci;i++){
    
    _chrmandlocus[i].resize(2);
    
    if (GetDistance(i) >= 100){
      cnum++;
      lnum=0;
      cstart[cnum] = i;
    } else if(cnum==-1){
      cerr << "first locus should have distance of >=100, but doesn't" << endl;
    }else{
      lnum++;
    }
    _chrmandlocus[i][0] = cnum;
    _chrmandlocus[i][1] = lnum;
    cfinish[cnum] = i;
  }
  
  Vector_d LengthOfChrm(cnum+1);
  NumberOfChromosomes = cnum +1;
  SizesOfChromosomes = new unsigned int[NumberOfChromosomes];
  
  Chromosome **C = new Chromosome*[cnum+1]; 
  for(int i=0;i<=cnum;i++){
    int size = cfinish[i] - cstart[i] + 1;
    C[i] = new Chromosome(size,cstart[i], populations);
    //C[i]->SetDistance(i,GetDistance(cstart[i]));
    stringstream label;
    string result;
    if( chrmlabels.size() == 0 ){
      label << "\"" << i << "\"";
      result = label.str();
      C[i]->SetLabel(0,result);
    }
    else
      C[i]->SetLabel(0,chrmlabels[cstart[i]]);
    for(int j=0;j<size;j++){
      (*(C[i]))(j) = (*this)(cstart[i]+j);
      C[i]->SetDistance(j,GetDistance(cstart[i]+j));
      if( j != 0 ){
	string s1("\"X\"");
	if( C[i]->GetLabel(0) != s1 ){
	  LengthOfGenome += GetDistance(cstart[i]+j);
	  LengthOfChrm(i) += GetDistance(cstart[i]+j);
	  //              cout << i << " " << j << " " << GetDistance(cstart[i]+j) << endl;
	}
	else{
	  LengthOfXchrm += GetDistance(cstart[i]+j);
	  X_data = true;
	  C[i]->ResetStuffForX();
	}
      }
    }
    SizesOfChromosomes[i] = C[i]->GetSize(); 
  }
  
  delete cstart;
  delete cfinish;
  return C;
}

CompositeLocus*& Genome::operator() ( int ElementNumber ) const
{
  if ( ElementNumber >= (int)NumberOfCompositeLoci ){
    cout << "WARNING: Genome::operator() Element Number";
    cout << ElementNumber << " > NumberOfCompositeLoci: " << NumberOfCompositeLoci << " accessed." << endl;
  }
  if ( ElementNumber < 0){
    cout << "WARNING: Genome::operator() ElementNumber";
    cout << ElementNumber << " < 0 accessed." << endl;
  }

  return TheArray[ElementNumber];
}

//creates array of CompositeLocus objects
void Genome::SetNumberOfCompositeLoci( int NumberOfElements )
{
  if(NumberOfCompositeLoci > 0){
    delete TheArray;
  }
  NumberOfCompositeLoci = NumberOfElements;
  TheArray = new CompositeLocus*[ NumberOfCompositeLoci ];
  Distances.SetNumberOfElements( NumberOfCompositeLoci );
}

void Genome::SetDistance( int locus, float distance )
{
  Distances( locus ) = distance;
}

void Genome::SetSizes(){
  TotalLoci = 0;

  for(unsigned int i = 0; i < NumberOfCompositeLoci; i++ ){
    TotalLoci += TheArray[i]->GetNumberOfLoci();
  }

}

//Accessors
unsigned int Genome::GetNumberOfCompositeLoci()
{
  return NumberOfCompositeLoci;
}

unsigned int Genome::GetNumberOfChromosomes(){
  return NumberOfChromosomes;
}

//returns total number of simple loci
unsigned int Genome::GetTotalNumberOfLoci(){
  return TotalLoci;
}
//returns int array of chromosome sizes
unsigned int *Genome::GetSizesOfChromosomes(){
  return SizesOfChromosomes;
}

Vector Genome::GetDistances()
{
  return( Distances );
}

float Genome::GetDistance( int locus )
{
  return( Distances( locus ) );
}

//returns total numbers of states accross all comp loci
int Genome::GetNumberOfStates()
{
  int ret = 0;
  for(unsigned int i = 0; i < NumberOfCompositeLoci; i++ ){
    ret += TheArray[i]->GetNumberOfStates();
  }
  return ret;
}

vector<int> Genome::GetChrmAndLocus( int j ){
  return _chrmandlocus[j];
}

vector<vector<int > > Genome::GetChrmAndLocus(){
  return _chrmandlocus;
}
bool Genome::isX_data()
{
   return X_data;
}

double Genome::GetLengthOfGenome()
{
   return LengthOfGenome;
}
double Genome::GetLengthOfXchrm()
{
   return LengthOfXchrm;
}









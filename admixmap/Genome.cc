#include <stdlib.h>
#include <sstream>
#include "Genome.h"
#include "Chromosome.h"

using namespace std;


Genome::Genome()
{
  NumberOfCompositeLoci = 0;
  NumberOfChromosomes = 0;
  TheArray = 0;
  LengthOfGenome = 0;
  LengthOfXchrm = 0;
  TotalLoci = 0;
}

//used by Chromosome
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
//   if(NumberOfCompositeLoci > 0)
//     for(unsigned int i = 0; i < NumberOfCompositeLoci; i++){
//       delete TheArray[i];
//     }
  delete[] TheArray;
  delete[] SizesOfChromosomes; 
}

//sets the labels of all composite loci in TheArray using SetLabel function in CompositeLocus
//all loci within each composite locus are assigned the same label from labels vector
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

void Genome::loadAlleleStatesAndDistances( vector<string> *ChrmLabels, AdmixOptions *options,InputData *data_, LogWriter *Log){
  string *LociLabelsCheck = 0;

  // Load number of allelic states and distances.
  //Also creates the CompositeLocus array
  Log->logmsg(false,"Loading ");
  Log->logmsg(false,options->getGeneInfoFilename());
  Log->logmsg(false,".\n");

  Matrix_d& locifileData = (Matrix_d&) data_->getGeneInfoMatrix();
  
  LociLabelsCheck = new string[ locifileData.GetNumberOfRows() ];

  //determine number of composite loci
  NumberOfCompositeLoci = locifileData.GetNumberOfRows();
    for( int i = 0; i < locifileData.GetNumberOfRows(); i++ )
     if( locifileData( i, 1 ) == 0.0 ) NumberOfCompositeLoci--;

  //set up CompositeLocus objects
  InitialiseCompositeLoci();

  //load locusfile data
  Log->logmsg(false,"Loading ");
  Log->logmsg(false,options->getGeneticDataFilename());
  Log->logmsg(false,".\n");

  // Set number of alleles at each locus
  int index =0;
  size_t next_line = 0;
  for(unsigned int i = 0; i < NumberOfCompositeLoci; i++ ){
    ++next_line;

    const Vector_s& m = data_->getGeneInfoData()[next_line];

    //set numbers of alleles and distances for each locus
    LociLabelsCheck[index] = m[0];
    if (m.size() == 4)
      ChrmLabels->push_back(m[3]);
    TheArray[i]->SetNumberOfAllelesOfLocus( 0, (int)locifileData( i, 0 ) );
    SetDistance( i, locifileData( index, 1 ) );
    while( index < locifileData.GetNumberOfRows() - 1 && locifileData( index + 1, 1 ) == 0 ){
      ++next_line;

      TheArray[i]->AddLocus( (int)locifileData( index + 1, 0 ) );
      index++;
      LociLabelsCheck[index] = m[0];
    }

    TheArray[i]->SetNumberOfLabels();
    index++;
    //Log->logmsg(false,(*Loci)(i)->GetNumberOfLoci());
    //Log->logmsg(false," ");
  }
  Log->logmsg(false,"\n");

  // checks of input data files should be in class InputData
  if( options->getTextIndicator() ){

    Vector_s labels = data_->getGeneticData()[0];

    Vector_d vtemp = locifileData.GetColumn(1);
    Log->logmsg(true, vtemp.GetNumberOfElements());Log->logmsg(true," simple loci\n");
    vtemp.AddElement(0); // Forces SetLabels method to ignore first row of loci.txt (GenotypesFile)
    // Add a sex column if it is not included
    if( ! options->genotypesSexColumn() ){
      labels.insert(labels.begin(), "\"extracol\"");
    }
    SetLabels(labels, vtemp);
  }

  index = 0;
  for(unsigned int i = 0; i < NumberOfCompositeLoci; i++ ){
    for( int j = 0; j < TheArray[i]->GetNumberOfLoci(); j++ ){
      if( LociLabelsCheck[index].compare( TheArray[i]->GetLabel(0) ) ){
	Log->logmsg(true, "Error in loci names in genotypes file and loci file at loci\n" );
	Log->logmsg(true, i);
	Log->logmsg(true, LociLabelsCheck[index] );
	Log->logmsg(true, " " );
	Log->logmsg(true, j);
	Log->logmsg(true, TheArray[i]->GetLabel(j) );
	Log->logmsg(true, "\n" );
	exit(0);
      }
      index++;
    }
  }  
  delete [] LociLabelsCheck;
}

//Creates an array of pointers to Chromosome objects, sets their labels
//also determines length of genome, NumberOfCompositeLoci, TotalLoci, NumberOfChromosomes, SizesOfChromosomes, 
//LengthOfXChrm and chrmandlocus
Chromosome** Genome::GetChromosomes( int populations, std::vector<std::string> &chrmlabels )
{
  int *cstart = new int[NumberOfCompositeLoci];
  int *cfinish = new int[NumberOfCompositeLoci];
  //since we don't know the number of chromsomes yet, these arrays have the maximum number, num. Comp. Loci, as size

  int cnum = -1; //number of chromosomes -1
  int lnum = 0;
  _chrmandlocus.resize(NumberOfCompositeLoci);
  X_data = false;
  
  TotalLoci = 0;
  for(unsigned int i = 0; i < NumberOfCompositeLoci; i++){
    
    _chrmandlocus[i].resize(2);
    
    if (GetDistance(i) >= 100){//new chromosome
      cnum++;
      lnum = 0; 
      cstart[cnum] = i; //locus number of first locus on new chromosome
    } else if(cnum==-1){
      cerr << "first locus should have distance of >=100, but doesn't" << endl;
    }else{
      lnum++;
    }
    _chrmandlocus[i][0] = cnum;
    _chrmandlocus[i][1] = lnum;
    cfinish[cnum] = i;//locus number of last locus on currrent chromosome
  }
  
  double LengthOfChrm[cnum+1];
  NumberOfChromosomes = cnum +1;
  SizesOfChromosomes = new unsigned int[NumberOfChromosomes];//array to store lengths of the chromosomes
  
  Chromosome **C = new Chromosome*[cnum+1]; 
  //C is an array of chromosome pointers

  for(int i = 0; i <= cnum; i++){//loop over chromsomes
    int size = cfinish[i] - cstart[i] + 1;//number of loci on chromosome i
    C[i] = new Chromosome(size, cstart[i], populations);
    //C[i] is a pointer to Chromosome

    //set chromosome label
    stringstream label;
    string result;
    //default (if none supplied)
    if( chrmlabels.size() == 0 ){
      label << "\"" << i << "\"";
      result = label.str();
      C[i]->SetLabel(result);
    }
    else
      C[i]->SetLabel(chrmlabels[cstart[i]]);

    for(int j = 0; j < size; j++){//loop over loci on chromosome
      //assign pointers to composite locus objects
      //NB: C[i] is a chromosome pointer, *(C[i]) is the chromosome it points to, 
      //(*(C[i]))(j) is an element of the array (of pointers) TheArray ie a locus on the chromosome
      //similarly, (*this)(j) is a locus on the genome
      //LHS and RHS are separate pointers but point to the SAME CompositeLocus objects
      //ie pointers are duplicated but CompositeLocus objects are not.

      (*(C[i]))(j) = (*this)(cstart[i]+j);
      C[i]->SetDistance(j,GetDistance(cstart[i]+j));

      if( j != 0 ){
	string s1("\"X\"");
	if( C[i]->GetLabel(0) != s1 ){
	  LengthOfGenome += GetDistance(cstart[i]+j);
	  LengthOfChrm[i] += GetDistance(cstart[i]+j);
	  //              cout << i << " " << j << " " << GetDistance(cstart[i]+j) << endl;
	  //NB length of genome and LengthOfChrm do not include X chromosome
	}
	//case of X chromosome
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

//accesses a composite locus
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
void Genome::InitialiseCompositeLoci()
{
  TheArray = new CompositeLocus*[ NumberOfCompositeLoci ];
  Distances.SetNumberOfElements( NumberOfCompositeLoci );
  for(unsigned int i = 0; i < NumberOfCompositeLoci;i++){
    //each element in the array is a pointer to a CompositeLocus object
    TheArray[i] = new CompositeLocus();
    TheArray[i]->SetNumberOfLoci(1);//not necessary
  }
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

int Genome::getNumberOfLoci(int j){
  return TheArray[j]->GetNumberOfLoci();
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









#include "Genome.h"
#include "Chromosome.h"
#include <stdlib.h>

using namespace std;

// composite-level methods

Genome::Genome()
{
  NumberOfCompositeLoci = 0;
  TheArray = 0;
  LengthOfGenome = 0;
  LengthOfXchrm = 0;
}

Genome::Genome( int NumberOfElements )
{
  if ( NumberOfElements < 1 ){
    NumberOfElements = 1;
  }
  NumberOfCompositeLoci = NumberOfElements;
  TheArray = new AbstractCompLocus*[ NumberOfCompositeLoci ];
  Distances.SetNumberOfElements( NumberOfCompositeLoci );
  LengthOfGenome = 0;
  LengthOfXchrm = 0;
}


Genome::~Genome()
{
  if(NumberOfCompositeLoci > 0){
    for(int i=0; i<NumberOfCompositeLoci; i++){
      //delete TheArray[i];
    }
    delete TheArray;
  }
}

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

void
Genome::accept(LocusVisitor& v)
{
  v.visitGenome(*this);
  for( int i = 0; i < size(); i++ ){
    (*this)(i)->accept(v);
  }
}

int
Genome::size()
{
  return NumberOfCompositeLoci;
}

vector<int> Genome::GetChrmAndLocus( int j ){
  return _chrmandlocus[j];
}

Genome*
Genome::GetChromosomes( int populations, vector<string> chrmlabels )
{
  int *cstart = new int[NumberOfCompositeLoci];
  int *cfinish = new int[NumberOfCompositeLoci];
  int cnum = -1;
  int lnum=0;
  _chrmandlocus.resize(NumberOfCompositeLoci);
  X_data = false;

  for(int i=0;i<NumberOfCompositeLoci;i++){
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
  Genome* arr = new Genome(cnum+1);
  for(int i=0;i<=cnum;i++){
     int size = cfinish[i] - cstart[i] + 1;
     Chromosome* chromosome = new Chromosome(size,cstart[i], populations);
     arr->SetDistance(i,GetDistance(cstart[i]));
     stringstream label;
     string result;
     if( chrmlabels.size() == 0 ){
        label << "\"" << i << "\"";
        result = label.str();
        chromosome->SetLabel(0,result);
     }
     else
        chromosome->SetLabel(0,chrmlabels[cstart[i]]);
     for(int j=0;j<size;j++){
        (*chromosome)(j) = (*this)(cstart[i]+j);
        chromosome->SetDistance(j,GetDistance(cstart[i]+j));
        if( j != 0 ){
           string s1("\"X\"");
           if( chromosome->GetLabel(0) != s1 ){
              LengthOfGenome += GetDistance(cstart[i]+j);
              LengthOfChrm(i) += GetDistance(cstart[i]+j);
//              cout << i << " " << j << " " << GetDistance(cstart[i]+j) << endl;
           }
           else{
              LengthOfXchrm += GetDistance(cstart[i]+j);
              X_data = true;
              chromosome->ResetStuffForX();
           }
        }
     }
    (*arr)(i) = chromosome;
  }
  
//  cout << LengthOfChrm << endl;
  delete cstart;
  delete cfinish;
  return arr;
}

bool
Genome::isX_data()
{
   return X_data;
}


double
Genome::GetLengthOfGenome()
{
   return LengthOfGenome;
}
double
Genome::GetLengthOfXchrm()
{
   return LengthOfXchrm;
}

AbstractCompLocus*&
Genome::operator() ( int ElementNumber ) const
{
  if ( ElementNumber >= NumberOfCompositeLoci ){
    cout << "WARNING: Genome::operator() Element Number";
    cout << ElementNumber << " > NumberOfCompositeLoci: " << NumberOfCompositeLoci << " accessed." << endl;
  }
  if ( ElementNumber < 0){
    cout << "WARNING: Genome::operator() ElementNumber";
    cout << ElementNumber << " < 0 accessed." << endl;
  }
  return TheArray[ ElementNumber ];
}

void
Genome::SetNumberOfCompositeLoci( int NumberOfElements )
{
  if(NumberOfCompositeLoci > 0){
    delete TheArray;
  }
  NumberOfCompositeLoci = NumberOfElements;
  TheArray = new AbstractCompLocus*[ NumberOfCompositeLoci ];
  Distances.SetNumberOfElements( NumberOfCompositeLoci );
}

int
Genome::GetNumberOfCompositeLoci()
{
  return NumberOfCompositeLoci;
}

void
Genome::SetDistance( int locus, float distance )
{
  Distances( locus ) = distance;
}

Vector
Genome::GetDistances()
{
  return( Distances );
}

float
Genome::GetDistance( int locus )
{
  return( Distances( locus ) );
}

// individual-level methods

void
Genome::AddLocus(int num)
{
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    TheArray[i]->AddLocus(num);
  }
}

Matrix_d
Genome::GetAlleleFreqs()
{
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    TheArray[i]->GetAlleleFreqs();
  }
  
  Matrix_d ret;
  return ret;
}

Vector_d
Genome::GetFst()
{
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    TheArray[i]->GetFst();
  }
  Vector_d ret;
  return ret;
}

Matrix_d
Genome::GetInfo()
{
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    TheArray[i]->GetInfo();
  }
  Matrix_d ret;
  return ret;
}

string
Genome::GetLabel(int index)
{
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    TheArray[i]->GetLabel(index);
  }
  string ret;
  return ret;
}

Vector_i
Genome::GetHapLabels(int index)
{
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    TheArray[i]->GetHapLabels(index);
  }
  Vector_i ret;
  return ret;
}

Matrix_i
Genome::GetLikelihoodAlleleFreqs()
{
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    TheArray[i]->GetLikelihoodAlleleFreqs();
  }
  Matrix_i ret;
  return ret;
}

int
Genome::GetMergedHaplotype(int num)
{
  int ret = 0;
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    ret += TheArray[i]->GetMergedHaplotype(num);
  }
  return ret;
}

int
Genome::GetNumberOfAllelesOfLocus(int num)
{
  int ret = 0;
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    ret += TheArray[i]->GetNumberOfAllelesOfLocus(num);
  }
  return ret;
}

int
Genome::GetNumberOfLoci()
{
  int ret = 0;
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    ret += TheArray[i]->GetNumberOfLoci();
  }
  return ret;
}

int
Genome::GetNumberOfMergedHaplotypes()
{
  int ret = 0;
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    ret += TheArray[i]->GetNumberOfMergedHaplotypes();
  }
  return ret;
}

int
Genome::GetNumberOfStates()
{
  int ret = 0;
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    ret += TheArray[i]->GetNumberOfStates();
  }
  return ret;
}

Vector_d
Genome::GetPriorAlleleFreqs(int num)
{
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    TheArray[i]->GetPriorAlleleFreqs(num);
  }
  Vector_d ret;
  return ret;
}

Matrix_d
Genome::GetScore()
{
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    TheArray[i]->GetScore();
  }
  Matrix_d ret;
  return ret;
}

Matrix_d
Genome::GetScoreSq()
{
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    TheArray[i]->GetScoreSq();
  }
  Matrix_d ret;
  return ret;
}

int
Genome::GetSize()
{
  int size = 0;
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    size += TheArray[i]->GetSize();
  }
  return size;
}

Vector_d
Genome::GetStatsForEta(int num)
{
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    TheArray[i]->GetStatsForEta(num);
  }
  Vector_d ret;
  return ret;
}

Matrix_d
Genome::GetSumAlleleFreqs()
{
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    TheArray[i]->GetSumAlleleFreqs();
  }
  Matrix_d ret;
  return ret;
}

int
Genome::IsRandom()
{
  int ret = 0;
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    ret += TheArray[i]->IsRandom();
  }
  return ret;
}

void
Genome::ResetLikelihoodAlleleFreqs()
{
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    TheArray[i]->ResetLikelihoodAlleleFreqs();
  }
}

void
Genome::ResetScoreForMisSpecOfAlleleFreqs()
{
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    TheArray[i]->ResetScoreForMisSpecOfAlleleFreqs();
  }
}

void
Genome::ResetSumAlleleFreqs()
{
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    TheArray[i]->ResetSumAlleleFreqs();
  }
}

void
Genome::SampleAlleleFreqs(int num)
{
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    TheArray[i]->SampleAlleleFreqs(num);
  }
}

void
Genome::SamplePriorAlleleFreqs(Vector_d eta)
{
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    TheArray[i]->SamplePriorAlleleFreqs(eta);
  }
}

void
Genome::SetAlleleFreqs(Matrix_d mat)
{
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    TheArray[i]->SetAlleleFreqs(mat);
  }
}

void
Genome::SetDefaultAlleleFreqs(int num)
{
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    TheArray[i]->SetDefaultAlleleFreqs(num);
  }
}

void
Genome::SetDefaultMergeHaplotypes(Vector_d alpha)
{
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    TheArray[i]->SetDefaultMergeHaplotypes(alpha);
  }
}

void
Genome::SetLabel(int index, string str)
{
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    TheArray[i]->SetLabel(index, str);
  }
}

void
Genome::SetHistoricalAlleleFreqs(Matrix_d mat)
{
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    TheArray[i]->SetHistoricalAlleleFreqs(mat);
  }
}

void
Genome::SetNumberOfAllelesOfLocus(int num1, int num2)
{
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    TheArray[i]->SetNumberOfAllelesOfLocus(num1,num2);
  }
}

void
Genome::SetNumberOfLoci(int num)
{
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    TheArray[i]->SetNumberOfLoci(num);
  }
}

void
Genome::SumScoreForMisSpecOfAlleleFreqs()
{
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    TheArray[i]->SumScoreForMisSpecOfAlleleFreqs();
  }
}

void
Genome::UpdateFst()
{
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    TheArray[i]->UpdateFst();
  }
}

void
Genome::UpdatePriorAlleleFreqsGlobal(int j, const vector<Vector_d>& mu)
{
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    TheArray[i]->UpdatePriorAlleleFreqs(j,mu[i]);
  }
}





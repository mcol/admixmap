#include "TuneRW.h"

TuneRW::TuneRW()
{
  k = 1;
}

TuneRW::TuneRW(int inw, double insigma0, double inmin, double inmax, double intarget)
{
  SetParameters(inw, insigma0, inmin, inmax, intarget);
}

TuneRW::~TuneRW()
{
}

void TuneRW::SetParameters(int inw, double insigma0, double inmin, double inmax, double intarget)
{
  w = inw;
  sigma0 = insigma0;
  sigma = sigma0;
  min = inmin;
  max = inmax;
  target = intarget;
  k = 1;
  count = 0;
  NumberAccepted = 0;
}

double TuneRW::GetSigma()
{
  return sigma;
}

double TuneRW::UpdateSigma(int NumberAccepted)
{
  double ProportionAccepted = (double)NumberAccepted / w;
  sigma = sigma + sigma0 * ( ProportionAccepted - target ) / k;
  if( sigma > max )
    sigma = max;
  else if( sigma < min )
    sigma = min;
  k++;
  return sigma;
}

void TuneRW::Event(bool accept)
{
  count++;
  NumberAccepted+=accept;
  if( count == w ){
    double ProportionAccepted = (double)NumberAccepted / w;
    sigma = sigma + sigma0 * ( ProportionAccepted - target ) / k;
    if( sigma > max )
      sigma = max;
    else if( sigma < min )
      sigma = min;
    k++;
    NumberAccepted = 0;
    count = 0;
  }  
}

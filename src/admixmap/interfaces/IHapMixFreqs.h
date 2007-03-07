#ifndef IHAPMIXFREQS_H_
#define IHAPMIXFREQS_H_

class IFreqArray;
class HapMixOptions;
class InputData;
class IGenome;
class LogWriter;

class IHapMixFreqs
{
public:
	// IHapMixFreqs();
	virtual ~IHapMixFreqs();
  virtual const IFreqArray& getHaploidGenotypeProbs() const = 0;
  virtual const IFreqArray& getDiploidGenotypeProbs() const = 0;
  virtual void Initialise(
      HapMixOptions* const options,
      InputData* const Data,
      IGenome *pLoci,
      LogWriter &Log) = 0;
};

#endif /*IHAPMIXFREQS_H_*/

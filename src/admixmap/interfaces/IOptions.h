#ifndef IOPTIONS_H_
#define IOPTIONS_H_

/**
 * Abstract class defining base Options class interface.
 * This class has two implementations:
 * Options for running the program
 * MockOptions for testing
 */

class IOptions
{
public:
  // virtual IOptions() = 0;
  virtual ~IOptions() = 0;
  virtual unsigned int getgenotypesSexColumn() const = 0;
  virtual bool getHWTestIndicator() const = 0;
  virtual bool isRandomMatingModel() const = 0;
  virtual int getPopulations() const = 0;
  virtual bool getHapMixModelIndicator() const = 0;
  virtual bool getTestForAllelicAssociation() const = 0;
};

#endif /*IOPTIONS_H_*/

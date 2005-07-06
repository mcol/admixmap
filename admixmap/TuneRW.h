// *-*-C++-*-*
#ifndef TuneRW_H
#define TuneRW_H 1

class TuneRW
{
public:
    TuneRW(int w, double sigma0, double min, double max, double target);
    TuneRW();
    ~TuneRW();

    void SetParameters(int w, double sigma0, double min, double max, double target);
    
    double GetSigma();
    double UpdateSigma(int NumberAccepted);
    void   Event(bool);

private:
   double sigma0; // Initial value of parameter of random walk being tuned.
   double sigma; // Current value of parameter of random walk being tuned.
   double min; // Minimum value of sigma
   double max; // Maximum value of sigma
   double target; // Target acceptance probability
   int k; // Number of times sigma sigma updated
   int w; // Frequency sigma is updated - 10 is good
   int count; // Number of iterations since sigma lsat updated
   int NumberAccepted; // Number of accepted proposals in random-walk since sgma last updated
};

#endif /* ! TuneRW_H */

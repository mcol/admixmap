// *-*-C++-*-*
#ifndef TuneRW_H
#define TuneRW_H 1

class TuneRW
{
public:
    TuneRW();
    TuneRW(int, double, double, double, double);
    ~TuneRW();

    void SetParameters(int, double, double, double, double);
    
    double GetSigma();
    double UpdateSigma(int);
    void   Event(bool);

private:
    double sigma0;
    double sigma;
    double min;
    double max;
    double target;
    int k;
    int w;
    int count;
    int NumberAccepted;
};

#endif /* ! TuneRW_H */

// *-*-C++-*-*
#ifndef ADMIXMAP_H
#define ADMIXMAP_H 1

#include <stdlib.h>    /* for exit, strtol */

#include <string.h>    /* for strcmp, strcpy */
#include <string>
#include "common.h"
#include "Latent.h"
#include "Regression.h"
#include "AlleleFreqs.h"
#include "AdmixOptions.h"
#include "StratificationTest.h"
#include "DispersionTest.h"
#include "ScoreTests.h"
#include "InputData.h"

#define ADMIXMAP_VERSION "1.6.1"

int main( int argc , char** argv );



#endif /* ADMIXMAP_H */

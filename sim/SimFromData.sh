 #!/bin/bash

export PRIORFILE=/opt/shared/admix/sarcoidosis/xchr/amass_gwammPrior.txt
export LOCUSFILE=/opt/shared/admix/sarcoidosis/xchr/amass_gwammLocus.txt
export GENOFILE=/opt/shared/admix/sarcoidosis/xchr/amass_gwammGeno.txt
export OUTCOMEFILE=/opt/shared/admix/sarcoidosis/xchr/amass_gwammOutcome.txt

R CMD BATCH --vanilla SimFromData.R SimFromDataR.log

admixmap simargs.txt
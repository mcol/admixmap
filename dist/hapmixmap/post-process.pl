#!/usr/bin/perl -w
use strict;

system("R CMD BATCH --quiet --no-save --no-restore ./AdmixmapOutput.R results/Rlog.txt");

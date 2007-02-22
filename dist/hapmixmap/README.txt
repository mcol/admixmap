HAPMIXMAP running instructions.

Prerequisites:
 * R, available from http://www.r-project.org/
Useful, but not essential:
 * Perl. Windows version available from http://www.activestate.com/Products/ActivePerl/

Once you have the prerequisites installed:

1. Unzip the hapmixmap file into some directory, for example C:\hapmixmap
2. Go to Start -> Run
3. Type "cmd" (without quotes) and press ENTER.
4. Type "cd hapmixmap" (without quotes) and press ENTER.

You are now ready to run the program. Use the perl script provided, which will run
the analysis for you.  The analysis consists of 2 steps:
1. Training with HapMap data only
2. Run with both HapMap data and case-control data and evaluate score test.

Each step can be resumed, to extend the number of iterations, after it is completed.
After each step, an R script is run to process the output.

To do step 1, type:                perl hapmixmap.pl --train.
To resume a training run, type: perl hapmixmap.pl --train --resume
To do step 2, type:                perl hapmixmap.pl --test
To resume a test run, type:    perl hapmixmap.pl --test --resume

You can also combine the two steps into one by typing: perl hapmixmap.pl --train --test.

If you do not have perl or prefer to run the program from the command-line, you can do so as follows:
To do step 1, type:                  ./hapmixmap training-initial.conf
To resume a training run, type: ./hapmixmap training-resume.conf
To do step2, type:                  ./hapmixmap testing.conf

If you do this, you will have to run the R script afterwards. First set the RESULTSDIR environment variable:
In Windows, type: set RESULTSDIR=Results
In Linux/UNIX, type: export RESULTSDIR=Results

Then to run the script:
R CMD BATCH --vanilla AdmixmapOutput.R Rlog.txt


HAPMIXMAP running instructions.

Prerequisites:

 * Perl, available from http://www.activestate.com/Products/ActivePerl/
 * R-project, available from http://www.r-project.com/

Once you have the prerequisites installed:

1. Unzip the hapmixmap file into some directory, for example C:\hapmixmap
2. Go to Start -> Run
3. Type "cmd" (without quotes) and press ENTER.
4. Type "cd hapmixmap" (without quotes) and press ENTER.

You are now ready to run the program. Type "run", like this:

C:\hapmixmap> run

This command will execute the provided "run.bat" script, which will run
the analysis for you.  The analysis consists of 4 steps. Each step will
be announced and will wait for you to press a key.

Each step can be executed separately, by typing a command.

I. Initial training
C:\hapmixmap> hapmixmap.exe training-initial.conf

II. Resume the training
C:\hapmixmap> hapmixmap.exe training-resume.conf

III. Test using the case-control data
C:\hapmixmap> hapmixmap.exe testing.conf

IV. Run the post-processing script
C:\hapmixmap> perl post-process.pl


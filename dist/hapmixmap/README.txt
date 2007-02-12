HAPMIXMAP running instructions.

First, open a DOS window:

1. Unzip the package into some directory, for example c:\hapmixmap
2. Go to Start -> Run
3. Type "cmd" (without quotes) and press ENTER.
4. Type "cd hapmixmap" (without quotes) and press ENTER.

You are now ready to run the program.

The analysis consists of 4 steps. Press ENTER after typing each line.


I. Initial training
> hapmixmap.exe training-initial.conf

II. Resume the training
> hapmixmap.exe training-resume.conf

III. Test using the case-control data
> hapmixmap.exe testing.conf

IV. Run the post-processing script
> perl post-process.pl


If you don't want to type all that by hand, there is a batch script
provided which runs all the commands for you. Just type "run" (without
quotes) and press ENTER to run it.

> run

The script will pause before each step, waiting for you to press ENTER.

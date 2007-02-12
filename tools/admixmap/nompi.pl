## template script for running independent parallel jobs without MPI. Provide rank and size as args
## NB no communication is possible inside the script

use strict;
use vars;

if($#ARGV >1 ){
    print("Too many arguments\n");
} 

if($#ARGV != 1){
    print("Please specify rank and number of processes\n\n");
    exit;
}

my $rank = $ARGV[0];
my $size = $ARGV[1];
if($rank >= $size){
    print("Invalid rank and size\n");
    exit;
}

print("I am $rank of $size\n");

#do something


#gather results
#? not wise as we don't know when each process is finished
my $reportfile = ">report.txt";
if($rank ==0){
    open (REPORTFILE, $reportfile) or die("Could not open report file\n");
}

if($rank == 0){
 close(REPORTFILE);   
}

print REPORTFILE "$rank\n";








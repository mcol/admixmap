rem This runs hapmixmap with tutorial configuration files, with the
rem diabetic nephropathy data set, containing 17 typed loci.

hapmixmap training-initial.conf
hapmixmap training-resume.conf
hapmixmap -ftesting.conf --ccgenotypesfile=data\CaseControlGenotypes_2_genotypes.txt --outcomevarfile=data\CaseControlGenotypes_2_outcome.txt

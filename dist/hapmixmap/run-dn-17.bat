rem This runs hapmixmap with tutorial configuration files, with the
rem diabetic nephropathy data set, containing 17 typed loci.

hapmixmap training-initial.conf
hapmixmap training-resume.conf
hapmixmap -ftesting.conf --ccgenotypesfile=data/CaseControlGenotypes2.txt --outcomevarfile=data/CaseControlOutcome2.txt

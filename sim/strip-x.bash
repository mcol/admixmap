#!/bin/bash


if ! [ -f data/genotypes.ped ]
    then
    echo "${0}: run SimAdmixture.R first.  Aborting." 1>&2
    exit 1
    fi


n_x=$(grep 'X$' data/loci.txt | wc -l)
n_not_x=$(grep -v 'X$' data/loci.txt | wc -l)


# Check for multiple executions
if true || [ "${n_x}" -ne 0 ]
    then

    sed -i "$((${n_not_x}+1)),\$d" data/loci.txt

    sed -i "$((${n_not_x}*2)),\$d" data/priorallelefreqs.txt

    mv data/genotypes.ped data/genotypes.ped.with-x

    gawk -v first_x=$((${n_not_x}+6)) \
	'
	# Every line
	    {
	    printf "%s", $1
	    for ( f=2 ; f<first_x ; ++f )
		printf " %s", $f
	    printf "\n"
	    }

	' >data/genotypes.ped <data/genotypes.ped.with-x

    fi

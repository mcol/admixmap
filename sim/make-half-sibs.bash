#!/bin/bash


if ! [ -f data/genotypes.ped ]
    then
    echo "${0}: run SimAdmixture.R first.  Aborting." 1>&2
    exit 1
    fi


head -6 data/genotypes.ped | \
    gawk '{ if (NR!=1) $1=1 ; if ($2==1||$2==2) $5=1 ; if ($2==3) $5=2 ; if ($2==4) {$3=1;$4=3} ; if ($2==5) {$3=2;$4=3} ; if ($2==1||$2==2||$2==3) for (f=7;f<NF;++f) $f="0,0" ; print $0 }' >| data/half-sibs.ped


echo "Topology summary of data/half-sibs.ped:"
gawk '{print $1 " " $2 " " $3 " " $4 " " $5 " " $6 " " $7}' < data/half-sibs.ped

#!/bin/bash
# Testing script

rm -rf bromba
perl start-project.pl \
	--name bromba \
	--population Eur \
	--genotypes ../../../hapmap/example-data/data/rawgenotypes.txt
tree bromba
cat bromba/project.ini

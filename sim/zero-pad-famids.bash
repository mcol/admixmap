#!/bin/bash

sed -i 's/^\([0-9]\)\t/F00\1\t/;s/^\([0-9][0-9]\)\t/F0\1\t/;s/^\([0-9][0-9][0-9]\)\t/F\1\t/' data/genotypes.txt data/genotypes.ped

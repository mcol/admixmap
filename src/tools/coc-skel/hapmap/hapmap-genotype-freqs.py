#!/usr/bin/python
# -*- coding: UTF-8 -*-
# 
# Generating R file with loci frequencies, based on a locus file in
# hapmap format and population name. Used as prior in the information
# reward calculation.
#
# Copyright (C) 2007 Maciej Bliziński
# 
# This script is free software distributed WITHOUT ANY WARRANTY. 
# You can redistribute it and/or modify it under the terms of the
# GNU General Public License, 
# version 2 or later, as published by the Free Software Foundation. 
# See the file COPYING for details.
# 
 
# Typical usage:
# python hapmap-genotype-freqs.py \
# --loci /home/maciej/chr22-data/Eur/chr22data/mi_loci.txt -p ceu \
# > ceu-prior.R

import psycopg

from userdata import UserData, RArray, LociPriors
from userdata import HapmixmapMaskingIndex
from optparse import OptionParser

def main():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option("-l", "--loci", dest = "locus_fn",
            help="Locus file", metavar = "FILE")
    parser.add_option("-p", "--population", dest = "population",
            help="Population name: ceu, yri, jpt-chb", metavar = "POP")
    parser.add_option("-n", "--index", dest = "in_index",
            help="Index file with the masked loci, mi_cc_index.txt", metavar = "FILE")

    (options, args) = parser.parse_args()
    assert options.locus_fn, "Locus file name needed"
    assert options.population, "Population name needed"
    assert options.in_index, "hapmixmap index file (mi_cc_index.txt) needed"

    mask = HapmixmapMaskingIndex(options.in_index)
    ud = UserData()
    ud.read_loci(options.locus_fn)

    conn = psycopg.connect("dbname=genepi")
    lp = LociPriors(mask.filter_loci(ud.loci), options.population, conn)
    print lp.format_dget()

if __name__ == '__main__':
    main()

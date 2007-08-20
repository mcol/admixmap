#!/usr/bin/python
# -*- coding: UTF-8 -*-
#
# Preparing data for IMPUTE
# http://www.stats.ox.ac.uk/%7Emarchini/software/gwas/impute.html
#
# This file is part of Genepi, genetic data analysis software.
# Copyright (C) 2007 Maciej Bliziński
# 
# Genepi is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# Genepi is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with Genepi; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

import sys
from optparse import OptionParser
from userdata import *
import os.path
import psycopg

def main():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option("-d", "--data-dir", dest = "data_dir",
            help="Data directory", metavar = "DIR",
            default = "/home/maciej/chr22-data/Eur/chr22data")
    parser.add_option("-t", "--training", dest = "training_fn",
            help="Input training haplotypes file", metavar = "FILE",
            default = "mi_genotypes.txt")
    parser.add_option("-m", "--masked", dest = "masked_fn",
            help="Input testing diplotypes file", metavar = "FILE",
            default = "mi_cc_masked.txt")
    parser.add_option("-l", "--train-loci", dest = "loci_fn",
            help="Input locus file", metavar = "FILE",
            default = "mi_loci.txt")
    parser.add_option("-a", "--haplo", dest = "out_haplo_fn",
            help="Output haplo.txt file", metavar = "FILE",
            default = "haplo.txt")
    parser.add_option("-b", "--geno", dest = "out_geno_fn",
            help="Output geno.txt file", metavar = "FILE",
            default = "geno.txt")
    parser.add_option("-c", "--legend", dest = "legend_fn",
            help="Output legend.txt file", metavar = "FILE",
            default = "legend.txt")

    (options, args) = parser.parse_args()

    in_train_haplo_fn = os.path.join(options.data_dir, options.training_fn)
    in_train_loci_fn = os.path.join(options.data_dir, options.loci_fn)
    in_test_fn = os.path.join(options.data_dir, options.masked_fn)

    ud = UserData()
    ud.read_loci(in_train_loci_fn)
    ud.read_haplotypes(in_train_haplo_fn)

    tud = DiploidUserData()
    tud.read_loci(in_train_loci_fn)
    tud.read_diplotypes(open(in_test_fn, "r"))

    # Genotypes coding is being read from a database. Visit
    # http://genedc.sf.net/ to get the data.
    conn = psycopg.connect("dbname='genepi'")

    imf = ImputeFormatter()
    imf.write_haplo(open(options.out_haplo_fn, "w"), ud.individuals, ud.loci)
    imf.write_geno(open(options.out_geno_fn, "w"), tud.individuals, tud.loci, conn)
    imf.write_legend(open(options.legend_fn, "w"), tud.loci, conn)

    # mi_genotypes.txt -> haplo.txt
    # mi_cc_masked.txt + mi_loci.txt -> geno.txt

if __name__ == '__main__':
    main()


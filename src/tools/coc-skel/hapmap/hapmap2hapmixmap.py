#!/usr/bin/python
# -*- coding: UTF-8 -*-
#
# Converts data from hapmap format to hapmixmap format.
#
# Copyright (C) 2007 Maciej Bliziński
# 
# This script is free software distributed WITHOUT ANY WARRANTY. 
# You can redistribute it and/or modify it under the terms of the
# GNU General Public License, 
# version 2 or later, as published by the Free Software Foundation. 
# See the file COPYING for details.
# 
 
import sys
from optparse import OptionParser
from userdata import *
import os.path

def main():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage)

    # Input files
    parser.add_option("-a", "--in-haplotypes", dest = "in_haplotypes",
            help="Input haplotypes file", metavar = "FILE")
    parser.add_option("-s", "--in-sample", dest = "in_sample",
            help="Input haplotypes file", metavar = "FILE")
    parser.add_option("-b", "--in-loci", dest = "in_loci",
            help="Input locus file", metavar = "FILE")
    parser.add_option("-c", "--out-haplotypes", dest = "out_haplotypes",
            help="Output haplotypes file", metavar = "FILE")
    parser.add_option("-d", "--out-loci", dest = "out_loci",
            help="Output locus file", metavar = "FILE")

    # Processing parameters
    parser.add_option("-o", "--output-dir", dest = "output_dir",
            help="Output directory, defaults to current.", metavar = "DIR")

    parser.set_defaults(output_dir = "")
    # Output files are hard-coded
    #
    # Files to be written:
    # http://actin.ucd.ie/trac/genepi/wiki/Coefficient%20of%20Constraint

    (options, args) = parser.parse_args()
    try:
        assert options.in_haplotypes, "Input haplotypes not specified"
        assert options.in_loci, "Input loci not specified"
        assert options.in_sample, "Sample not specified"
        assert options.out_haplotypes, "Output haplotypes not specified"
        assert options.out_loci, "Output loci not specified"
    except AssertionError, e:
        print parser.format_help()
        print
        print "Error:", e
        print
        sys.exit(1)

    ud = UserData()
    ud.read_hapmap_legend(options.in_loci)
    ud.read_hapmap_haplotypes(options.in_haplotypes, options.in_sample)
    print "Calculating loci to write"
    loci_to_write = filter(
            lambda x: not x.is_monomorphic(ud.individuals),
            ud.loci)

    # ud.randomize_masking(0, 0, 0)

    dmm = DataMaskingMediator(options.output_dir)
    hmf = HapmixmapFormatter()

    print "Writing haplotypes."
    out_haplo_fh = open(dmm.get_output_filename(options.out_haplotypes), "w")
    out_haplo_fh.write(hmf.format_header(
        ud.individuals, loci_to_write, []))
    out_haplo_fh.write(hmf.format_haplotypes(
        ud.individuals, loci_to_write, [], mask = False))
    out_haplo_fh.close()

    print "Writing loci."
    lff = LocusFileFormatter()
    mi_loci = open(dmm.get_output_filename(options.out_loci), "w")
    mi_loci.write(lff.format(loci_to_write))
    mi_loci.close()

if __name__ == '__main__':
    # print "This program isn't finished."
    # sys.exit(1)
    main()

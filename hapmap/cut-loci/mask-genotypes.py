#!/usr/bin/python
# -*- coding: UTF-8 -*-
#
# Rewrite of genotypes masking.
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

def main():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage)

    # Input files
    parser.add_option("-a", "--in-haplotypes", dest = "in_haplotypes",
            help="Input haplotypes file", metavar = "FILE")
    parser.add_option("-b", "--in-loci", dest = "in_loci",
            help="Input locus file", metavar = "FILE")

    # Processing parameters
    parser.add_option("-i", "--percent-individuals", dest = "percent_individuals", type = "int",
            help="Percent individuals to mask.")
    parser.add_option("-l", "--percent-loci", dest = "percent_loci", type = "int",
            help="Percent loci to mask.")
    parser.add_option("-m", "--limit-loci", dest = "limit_loci", type = "int",
            help="Limit results to given number of loci.")

    parser.add_option("-d", "--output-dir", dest = "output_dir",
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
        assert options.percent_individuals, "Percent individuals not specified"
        assert options.percent_loci, "Percent loci not specified"
    except AssertionError, e:
        print parser.format_help()
        print
        print "Error:", e
        print
        sys.exit(1)

    ud = UserData()
    ud.read_loci(options.in_loci)
    ud.read_haplotypes(options.in_haplotypes)
    ud.randomize_masking(options.percent_individuals,
            options.percent_loci, options.limit_loci)

    dmm = DataMaskingMediator(options.output_dir)
    #  mi_cc_fastphase.inp   with question marks     5000 loci   CC,
    #  masked    fastPHASE   ?
    fpf = FastPhaseFormatter()
    fpf.write(dmm.get_output_filename("mi_cc_fastphase.inp"), ud.get_masked_individuals(),
            ud.loci_to_write, ud.masked_loci)

    #  mi_cc_index.txt   maskedindivs = <vals> \n maskedloci = <vals>
    #  -     -   hapmixmap   -
    hif = HapmixmapIndexFormatter()
    hif.write(dmm.get_output_filename("mi_cc_index.txt"), ud.get_masked_individuals(),
            ud.loci_to_write, ud.masked_loci, ud.get_unmasked_individuals())

    #  mi_cc_masked.txt 5000 loci,   CC, masked      hapmixmap D
    hmf = HapmixmapFormatter()
    mi_cc_masked = open(dmm.get_output_filename("mi_cc_masked.txt"), "w")
    mi_cc_masked.write(hmf.format_header(
        ud.get_masked_individuals(), ud.loci_to_write, ud.masked_loci))
    mi_cc_masked.write(hmf.format_diplotypes(
        ud.get_masked_individuals(), ud.loci_to_write, ud.masked_loci, mask = True))
    mi_cc_masked.close()
    
    #  mi_cc_observed_dput.txt 1500 loci   CC, unmasked    R   D
    dgf = DgetFormatter()
    mi_cc_observed_dput = open(dmm.get_output_filename("mi_cc_observed_dput.txt"), "w")
    mi_cc_observed_dput.write(dgf.format(ud.get_masked_individuals(),
        ud.loci_to_write, ud.masked_loci))
    mi_cc_observed_dput.close()

    #  mi_cc_observed.txt 1500 loci   CC, unmasked    hapmixmap  D
    mi_cc_observed = open(dmm.get_output_filename("mi_cc_observed.txt"), "w")
    mi_cc_observed.write(hmf.format_header(
        ud.get_masked_individuals(), ud.loci_to_write, ud.masked_loci))
    mi_cc_observed.write(hmf.format_diplotypes(
        ud.get_masked_individuals(), ud.loci_to_write, ud.masked_loci, mask = False))
    mi_cc_observed.close()

    #  mi_cc_original_full.txt   the same as mi_cc.txt, except for the
    #  quotes around field names     5000 loci   CC, unmasked

    mi_cc_original_full = open(dmm.get_output_filename("mi_cc_original_full.txt"), "w")
    mi_cc_original_full.write(hmf.format_header(
        ud.get_masked_individuals(), ud.loci_to_write, ud.masked_loci))
    mi_cc_original_full.write(hmf.format_diplotypes(
        ud.get_masked_individuals(), ud.loci_to_write, ud.masked_loci, mask = False))
    mi_cc_original_full.close()

    #  hapmixmap     D
    #  mi_cc.txt         5000 loci   CC, unmasked    hapmixmap   D
    mi_cc = open(dmm.get_output_filename("mi_cc.txt"), "w")
    mi_cc.write(hmf.format_header(
        ud.get_masked_individuals(), ud.loci_to_write, ud.masked_loci))
    mi_cc.write(hmf.format_diplotypes(
        ud.get_masked_individuals(), ud.loci_to_write, ud.masked_loci, mask = False))
    mi_cc.close()

    #  mi_genotypes.txt          5000 loci   Training    hapmixmap   H
    mi_genotypes = open(dmm.get_output_filename("mi_genotypes.txt"), "w")
    mi_genotypes.write(hmf.format_header(
        ud.get_unmasked_individuals(), ud.loci_to_write, ud.masked_loci))
    mi_genotypes.write(hmf.format_haplotypes(
        ud.get_unmasked_individuals(), ud.loci_to_write, ud.masked_loci, mask = False))
    mi_genotypes.close()

    #  mi_loci.txt   locus file      5000 loci   -   hapmixmap   -
    lff = LocusFileFormatter()
    mi_loci = open(dmm.get_output_filename("mi_loci.txt"), "w")
    mi_loci.write(lff.format(ud.loci_to_write))
    mi_loci.close()

    #  mi_merged_cc_train.txt        5000 loci   Training + CC, masked
    #  hapmixmap     H + D
    mi_merged_cc_train = open(dmm.get_output_filename("mi_merged_cc_train.txt"), "w")
    mi_merged_cc_train.write(hmf.format_header(
        ud.get_unmasked_individuals(), ud.loci_to_write, ud.masked_loci))
    mi_merged_cc_train.write(hmf.format_haplotypes(
        ud.get_unmasked_individuals(), ud.loci_to_write, ud.masked_loci, mask = False))
    mi_merged_cc_train.write(hmf.format_diplotypes(
        ud.get_masked_individuals(), ud.loci_to_write, ud.masked_loci, mask = True))
    mi_merged_cc_train.close()

if __name__ == '__main__':
    main()

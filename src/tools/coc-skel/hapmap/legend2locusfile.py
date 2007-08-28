#!/usr/bin/python
# -*- coding: UTF-8 -*-
# 
# Hapmap.org to SQL converter
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
    parser.add_option("-a", "--hapmap", dest = "hapmap_dir",
            help="Input haplotypes file", metavar = "DIR")

    parser.set_defaults(output_dir = "")
    # /home/maciej/ucd/genedc-bzr/src/hapmap/download

    (options, args) = parser.parse_args()
    try:
        assert options.hapmap_dir, "Hapmap directory not specified"
    except AssertionError, e:
        print parser.format_help()
        print
        print "Error:", e
        print
        sys.exit(1)
    hapmap = Hapmap()
    for chr in hapmap.get_chromosomes():
        for pop in hapmap.get_populations():
            # Splitting data into chunks
            leg_fn = os.path.join(options.hapmap_dir,
                    hapmap.get_legend_fn(chr, pop))
            hap_fn = os.path.join(options.hapmap_dir,
                    hapmap.get_haplotypes_fn(chr, pop))
            sam_fn = os.path.join(options.hapmap_dir,
                    hapmap.get_sample_fn(chr, pop))
            ud = UserData()
            ud.read_hapmap_legend(os.path.join(os.path.join(
                options.hapmap_dir,
                "genotypes-%s-%s-r21-nr-fwd-legend.txt" % (chr, pop))))
            print "Writing locusfile for", chr, pop
            fh = open("%s-%s-loci.txt" % (chr, pop), "w")
            lff = LocusFileFormatter()
            fh.write(lff.format(ud.loci))
            fh.close()
            del lff
            del ud


if __name__ == '__main__':
    # import cProfile
    # cProfile.run("main()", "sqlprof")
    main()


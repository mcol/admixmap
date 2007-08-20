#!/usr/bin/python
# -*- coding: UTF-8 -*-
# 
# Hapmap.org to SQL converter
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
    import psycopg
    conn = psycopg.connect("dbname='genepi'")
    loci_numbers = {}
    for chr in hapmap.get_chromosomes():
    # for chr in ['chr22']:
        for pop in hapmap.get_populations():
            print pop
            leg_fn = os.path.join(options.hapmap_dir,
                    hapmap.get_legend_fn(chr, pop))
            hap_fn = os.path.join(options.hapmap_dir,
                    hapmap.get_haplotypes_fn(chr, pop))
            sam_fn = os.path.join(options.hapmap_dir,
                    hapmap.get_sample_fn(chr, pop))
            si = SqlInserter(conn, hap_fn, leg_fn, sam_fn,
                    loci_numbers)
            si.write(pop)

if __name__ == '__main__':
    # import cProfile
    # cProfile.run("main()", "sqlprof")
    main()


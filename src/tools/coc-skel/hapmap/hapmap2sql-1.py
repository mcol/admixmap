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
    import psycopg
    conn = psycopg.connect("dbname='genepi'")
    curs = conn.cursor()
    print "Building a dictionary of loci"
    curs.execute("SELECT id, snp_id FROM locus;")
    loci_numbers = {}
    for row in curs.fetchall():
        loci_numbers[row[1]] = row[0]
    curs.close()
    for chr in hapmap.get_chromosomes():
        for pop in hapmap.get_populations():
            # Splitting data into chunks
            leg_fn = os.path.join(options.hapmap_dir,
                    hapmap.get_legend_fn(chr, pop))
            hap_fn = os.path.join(options.hapmap_dir,
                    hapmap.get_haplotypes_fn(chr, pop))
            sam_fn = os.path.join(options.hapmap_dir,
                    hapmap.get_sample_fn(chr, pop))
            sc = SqlConverter(conn, hap_fn, leg_fn, sam_fn,
                    loci_numbers)
            print "Writing SQL for", chr, pop
            fh = open("%s-%s.sql" % (chr, pop), "w")
            sc.write(fh, pop)
            fh.close()
            del sc


if __name__ == '__main__':
    # import cProfile
    # cProfile.run("main()", "sqlprof")
    main()


#!/usr/bin/python
# -*- coding: UTF-8 -*-
# 
# A script to extract posterior probabilities from fastPHASE output.
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
import re

from optparse import OptionParser

class FastPhasePosteriorProbabilityDistribution:
    GENOTYPES = ("1,1", "1,2", "2,2")
    def __init__(self, fn, no_indivs):
        no_indivs = int(no_indivs)
        assert type(no_indivs) == int
        lines = open(fn, "r").readlines()
        lines.reverse()
        sample_re = re.compile(r"Sample from")
        self.individuals = {}
        sample_no = -1
        line_within_sample = 0
        while len(lines):
            line = lines.pop().strip()
            line_within_sample += 1;
            if re.search(sample_re, line):
                sample_no += 1
                line_within_sample = 0
                # print "Sample ", sample_no
                sys.stdout.write(".")
                sys.stdout.flush()
                continue
            if line_within_sample >= (no_indivs * 2):
                line_within_sample += 1
                continue
            # This must be an individual. Let's fetch his next line
            line2 = lines.pop().strip()
            line_within_sample += 1;
            indiv_no = line_within_sample / 2 - 1
            # print "Individual no", indiv_no
            if not self.individuals.has_key(indiv_no):
                self.individuals[indiv_no] = Individual(indiv_no, self)
            ind = self.individuals[indiv_no]
            # print "Updating the individuals, lines:"
            # print line[:72]
            # print line2[:72]
            ind.update(sample_no, (line, line2))
        print
    def test(self):
        print "Number of individuals:", len(self.individuals.keys())
        for ind in self.individuals.keys()[0:5]:
            self.individuals[ind].test()
    def dput(self, fn, masked_loci = None):
        if masked_loci:
            assert type(masked_loci) == list
        f = open(fn, "w")
        f.write("structure(.Data = c(\n")
        # Outermost: locus
        # Then: individual
        # Innermost: genotype
        no_indivs = len(self.individuals.keys())
        no_loci = len(self.individuals[0].loci)
        if masked_loci:
            loci_idx = masked_loci
        else:
            loci_idx = range(no_loci)
        # print "loci_idx:", loci_idx[:10]
        for i_key in self.individuals.keys():
            f.write("# Individual %s\n" % i_key)
            for locus_no in loci_idx:
                f.write("# Locus %s\n" % locus_no)
                gd = self.individuals[i_key].get_genotype_distrib(locus_no)
                assert sum([gd[i] for i in self.GENOTYPES]) - 1 < 1e-5
                for genotype in self.GENOTYPES:
                    f.write("%s, " % gd[genotype])
                f.write("\n")
        f.write("),\n")
        f.write(".Dim = c(%s, %s, %s),\n" % (3, len(loci_idx), no_indivs))
        f.write(".Dimnames = list(c(\"Genotype1\", \"Genotype2\", \"Genotype3\"), 1:%s, 1:%s))\n" % (len(loci_idx), no_indivs))
        f.close()

class Individual:
    def __init__(self, id, data_set):
        self.loci = {}
        self.id = id
        self.data_set = data_set
        self.ignored = False
    def update(self, sample_no, lines):
        assert type(lines) == tuple
        assert len(lines) == 2
        assert len(lines[0]) == len(lines[1])
        no_loci = len(lines[0])
        for i in range(no_loci):
            if not self.loci.has_key(i):
                self.loci[i] = {}
                for g in self.data_set.GENOTYPES:
                    self.loci[i][g] = 0
            locus = self.loci[i]
            g = self.get_genotype((lines[0][i], lines[1][i]))
            if g is None:
                continue
            assert g in self.data_set.GENOTYPES
            locus[g] += 1
            # locus.append((lines[0][i], lines[1][i]))
    def __str__(self):
        return "I: %s ..." % self.loci[0]
    def test(self):
        print "Individual", self.id
        for loc in self.loci.keys()[:5]:
            print self.loci[loc]
    def get_genotype(self, a):
        assert type(a) == tuple
        assert len(a) == 2
        if "?" in a:
            return None
        b = list(a)
        b.sort()
        g = "%s,%s" % (b[0], b[1])
        try:
            assert g in self.data_set.GENOTYPES
        except AssertionError, e:
            print "Invalid genotypes,", g, "is not one of", self.data_set.GENOTYPES
            raise
        # print "get_genotype(): incoming %s, outgoing %s" % (a, g)
        return g
    def get_genotype_counts(self, locus_no):
        return self.loci[locus_no]
    def get_genotype_distrib(self, locus_no):
        counts = self.get_genotype_counts(locus_no)
        total = sum([counts[i] for i in counts.keys()])
        # print counts, total
        reciprocal = 1.0 / float(total)
        for key in counts.keys():
            counts[key] *= reciprocal
        assert sum([counts[i] for i in counts.keys()]) - 1 < 1e-5
        # assert len(counts.keys()) == 3
        return counts

def main():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option("-s", "--sampled", dest = "file_name",
            help="Input file name, usually fastphase_sampledHgivG.txt", metavar = "FILE")
    parser.add_option("-o", "--output", dest = "out_file_name",
            help="Output file name", metavar = "FILE")
    parser.add_option("-i", "--individuals", dest = "no_indivs", type = "int",
            help="Number of individuals")
    parser.add_option("-f", "--args-file", dest = "args_file",
            help="File with masked loci index. This is usually a file with hapmixmap arguments.")

    (options, args) = parser.parse_args()
    if not options.file_name:
        print "Please specify the input file."
        sys.exit()
    if not options.out_file_name:
        print "Please specify the output file."
        sys.exit()
    if not options.no_indivs:
        print "Please specify the number of individuals."
        sys.exit()
    fp = FastPhasePosteriorProbabilityDistribution(options.file_name, options.no_indivs)
    print "Writing the", options.out_file_name, "file."
    masked_loci = None
    if options.args_file:
        for line in open(options.args_file, "r").readlines():
            if re.search(r'maskedloci', line):
                masked_loci_line = line
        # print masked_loci_line[:72]
        masked_loci = re.findall(r'\w+', masked_loci_line.split("=")[1])
        # print masked_loci[:10]
        masked_loci = map(int, masked_loci)
        # print masked_loci[:10]
        masked_loci = map(lambda(x): x - 1, masked_loci)
        # print masked_loci[:10]
    fp.dput(options.out_file_name, masked_loci)

if __name__ == '__main__':
    main()

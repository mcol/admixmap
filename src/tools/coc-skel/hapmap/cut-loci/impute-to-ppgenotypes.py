#!/usr/bin/python
# -*- coding: UTF-8 -*-
# 
# Converting data from IMPUTE output format to R object with posterior
# probabilities, such as the one produced by HAPMIXMAP
# (PPGenotypeProbs.txt).
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
import re
import itertools

# Python 2.3 compatibility
from sets import Set as set

class PosteriorProbabilitiesRArray(RArray):
    def __init__(self, data, dims, labels, snp_ids):
        RArray.__init__(self, data, dims, labels)
        self.snp_ids = snp_ids
    def format_data(self, loci):
        """The data format is:
        { 'rsXXX': [[0.1, 0.2, 0.7], [ ... ]], ... }
        """
        # print "PosteriorProbabilitiesRArray::format_data()"
        prefix = ""
        r = ""
        cnt = itertools.count()
        # Loop nesting:
        # 1. Individuals
        # 2. SNPs
        # 3. Genotypes
        for indiv_idx in xrange(0, len(loci[self.snp_ids[0]])):
            r += "\n# indiv %s\n" % indiv_idx
            for snp_id in self.snp_ids:
                for prob in loci[snp_id][indiv_idx]:
                    r += prefix
                    n = cnt.next()
                    if n % 3 == 0:
                        r += "\n# %s\n" % snp_id
                    r += "%s" % prob
                    if not prefix:
                        prefix = ", "
        return r

class PosteriorProbabilities(object):
    """Holds a data structure for posterior probabilities."""
    def __init__(self, mask):
        assert isinstance(mask, HapmixmapMaskingIndex)
        self.snp_ids = []
        self.loci = {}
        self.mask = mask
    def add_locus(self, snp_id, grouped_probs):
        assert len(grouped_probs) > 0
        # Check if all groups have 3 elements
        assert [3] * len(grouped_probs) == [len(x) for x in grouped_probs]
        self.snp_ids.append(snp_id)
        self.loci[snp_id] = grouped_probs

    def dput(self, fn):
        open(fn, "w").write(self.format_dput())

    def format_dput(self):
        # print self.snp_ids
        try:
            indivs = len(self.loci[self.snp_ids[0]])
        except IndexError:
            print self.snp_ids
            raise
        # print self.loci[self.snp_ids[0]]
        # print indivs, "indivs"
        # plain_probs = reduce(lambda x, y: x + y, [self.loci[x] for x in self.snp_ids])
        a = PosteriorProbabilitiesRArray(self.loci,
                [3, len(self.snp_ids), indivs],
                [
                    ['1', '2', '3'],
                    self.snp_ids, ["i%s" % x for x in xrange(indivs)], ],
                self.snp_ids)
        return a.format()
    def get_plain_probs(self):
        """Create an iterator over plain loci."""
        return self.PlainProbs(self.loci, self.snp_ids)

class ImputePosteriorProbabilities(PosteriorProbabilities):
    """Reads IMPUTE output format."""
    line_format = r'^(?P<junk>\w+) (?P<snp_id>rs[0-9]+) (?P<pos>[0-9\.]+) (?P<a1>[ACGT]) (?P<a2>[ACGT]) (?P<probs>[0-9\.\s]+)$'
    def __init__(self, fh, mask):
        """Read IMPUTE output format:

rs1000112 rs1000112 12166546 C G 0 0 0 0
        """
        PosteriorProbabilities.__init__(self, mask)
        assert isinstance(fh, file), "fh needs to be a file"
        # Masking index are 1-based
        for lineno in itertools.count(1):
            line = fh.readline()
            if not line.strip(): break
            if not lineno in self.mask.idx: continue
            parsed = re.match(self.line_format, line)
            assert parsed, "Parsing error: '%s'" % line
            grps = parsed.groupdict()
            probs_list = re.split(r'\s+', grps['probs'].strip())
            if len(probs_list) % 3 != 0:
                 print "Number of probabilities must be divisible by 3"
                 print len(probs_list)
                 print probs_list
            # grouped_probs = map(lambda x: probs_list[x:x + 3],
            #         itertools.islice(range(len(probs_list)), 0, len(probs_list), 3))
            grouped_probs = group_list(probs_list, 3)
            self.add_locus(grps['snp_id'], grouped_probs)

def main():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option("-i", "--input", dest = "in_impute",
            help="Input file, from IMPUTE, usually named out.txt", metavar = "FILE")
    parser.add_option("-o", "--output", dest = "dput_fn",
            help="Output PPGenotypeProbs.txt file", metavar = "FILE",
            default = "PPGenotypeProbs.txt")
    parser.add_option("-n", "--index", dest = "in_index",
            help="Index file with the masked loci, mi_cc_index.txt", metavar = "FILE")
    (options, args) = parser.parse_args()

    assert options.in_impute, "IMPUTE file name needed"
    assert options.in_index, "hapmixmap index file (mi_cc_index.txt) needed"

    mask = HapmixmapMaskingIndex(options.in_index)
    pp = ImputePosteriorProbabilities(open(options.in_impute, "r"), mask)
    # print pp.loci
    # a = RArray([0, 1, 2, 3], [2, 2], [['a', 'b'], ['c', 'd']])
    open(options.dput_fn, "w").write(pp.format_dput())

if __name__ == '__main__':
    main()


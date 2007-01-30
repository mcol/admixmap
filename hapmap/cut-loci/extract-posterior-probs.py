#!/usr/bin/python

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
                print "New sample:", sample_no
                continue
            if line_within_sample >= (no_indivs * 2):
                line_within_sample += 1
                continue
            # It must be an individual. Let's fetch his next line
            line2 = lines.pop().strip()
            line_within_sample += 1;
            indiv_no = line_within_sample / 2 - 1
            print "Individual no", indiv_no
            if not self.individuals.has_key(indiv_no):
                self.individuals[indiv_no] = Individual(indiv_no)
            ind = self.individuals[indiv_no]
            # print "Updating the individuals, lines:"
            # print line[:72]
            # print line2[:72]
            ind.update(sample_no, (line, line2))
    def test(self):
        print "Number of individuals:", len(self.individuals.keys())
        for ind in self.individuals.keys()[0:5]:
            self.individuals[ind].test()
    def dput(self, fn):
        f = open(fn, "w")
        f.write("structure(.Data = c(\n")
        # Outermost: locus
        # Then: individual
        # Innermost: genotype
        no_indivs = len(self.individuals.keys())
        no_loci = len(self.individuals[0].loci)
        for locus_no in range(no_loci):
            f.write("# Locus %s\n" % locus_no)
            for i_key in self.individuals.keys():
                f.write("# Individual %s\n" % i_key)
                gd = self.individuals[i_key].get_genotype_distrib(locus_no)
                for genotype in self.GENOTYPES:
                    if gd.has_key(genotype):
                        val = gd[genotype]
                    else:
                        val = 0.0
                    f.write("%s, " % val)
                f.write("\n")
        f.write("),\n")
        f.write(".Dim = c(%s, %s, %s),\n" % (3, no_indivs, no_loci))
        f.write(".Dimnames = list(c(\"Genotype1\", \"Genotype2\", \"Genotype3\"), 1:%s, 1:%s))\n" % (no_indivs, no_loci))
        f.close()

class Individual:
    def __init__(self, id):
        self.loci = {}
        self.id = id
    def update(self, sample_no, lines):
        assert type(lines) == tuple
        assert len(lines) == 2
        assert len(lines[0]) == len(lines[1])
        no_loci = len(lines[0])
        for i in range(no_loci):
            if not self.loci.has_key(i):
                self.loci[i] = []
            locus = self.loci[i]
            locus.append((lines[0][i], lines[1][i]))
    def __str__(self):
        return "I: %s ..." % self.loci[0]
    def test(self):
        print "Individual", self.id
        for loc in self.loci.keys()[:5]:
            print self.loci[loc]
    def get_genotype(self, a):
        assert type(a) == tuple
        assert len(a) == 2
        b = list(a)
        b.sort()
        g = "%s,%s" % (b[0], b[1])
        # print "get_genotype(): incoming %s, outgoing %s" % (a, g)
        return g
    def get_genotype_counts(self, locus_no):
        locus = self.loci[locus_no]
        counts = {}
        for sample in locus:
            g = self.get_genotype(sample)
            if not counts.has_key(g):
                counts[g] = 1
            else:
                counts[g] += 1
        return counts
    def get_genotype_distrib(self, locus_no):
        counts = self.get_genotype_counts(locus_no)
        total = sum([counts[i] for i in counts.keys()])
        # print counts, total
        for key in counts.keys():
            counts[key] /= float(total)
        return counts

def main():
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    parser.add_option("-s", "--sampled", dest = "file_name",
            help="Input file name, usually fastphase_sampledHgivG.txt", metavar = "FILE")
    parser.add_option("-o", "--output", dest = "out_file_name",
            help="Output file name", metavar = "FILE")
    parser.add_option("-i", "--individuals", dest = "no_indivs",
            help="Number of individuals")

    (options, args) = parser.parse_args()
    # print options, args
    # print "in: ", options.file_name
    # print "out: ", options.out_file_name
    # print "indivs: ", options.no_indivs
    fp = FastPhasePosteriorProbabilityDistribution(options.file_name, options.no_indivs)
    fp.dput(options.out_file_name)

if __name__ == '__main__':
    main()

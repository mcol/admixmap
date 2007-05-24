#!/usr/bin/python
# -*- coding: UTF-8 -*-
#
# Python classes handling genetic data
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

import re
import csv
import sys
import random
import os.path

class Locus:
    def __init__(self, snp_id, num_alleles, distance, previous_position, idx):
        self.idx = idx
        self.masked_idx = None
        self.snp_id = snp_id
        self.num_alleles = num_alleles
        self.distance = distance
        if distance:
            self.position = long(previous_position + 1e6 * distance)
        else:
            self.position = 0
    def is_monomorphic(self, individuals):
        assert type(individuals) == list
        last_genotype = individuals[0].get_genotype(self.snp_id)
        flipped = False
        for indiv in individuals[1:]:
            current_genotype = indiv.get_genotype(self.snp_id)
            if last_genotype != current_genotype:
                flipped = True
                break
            last_genotype = current_genotype
        return not flipped
    def set_masked_idx(self, idx):
        # print "Setting masked_idx of locus", self.snp_id, "(idx is", self.idx, ") to", idx
        self.masked_idx = idx

class Individual:
    def __init__(self, id, idx):
        self.idx = idx
        self.genotypes = {}
        self.id = id
        # print "Individual(%s)" % repr(id)
    def add_genotype(self, locus_id, genotype):
        # print "add_genotype(%s, %s)" % (locus_id, genotype)
        self.genotypes[locus_id] = genotype
    def getId(self):
        return self.id
    def get_genotype(self, snp_id):
        if self.genotypes.has_key(snp_id):
            return self.genotypes[snp_id]
        else:
            # Perhaps throw an exception?
            print "Snp %s not found in individual %s." % \
                    (snp_id, self.id)
            sys.exit(1)

class UserData:
    def __init__(self):
        self.individuals = []
        self.individuals_by_id = {}
        self.loci = []
        self.loci_by_id = {}
        self.whitespace = re.compile("\s+")
        self.masked_loci = {}
        self.masked_indivs = {}
        self.loci_to_write = []
    def get_unmasked_individuals(self):
        return filter(lambda x: not self.masked_indivs.has_key(x.id),
                self.individuals)
    def get_masked_individuals(self):
        return filter(lambda x: self.masked_indivs.has_key(x.id),
                self.individuals)
    def add_individual(self, individual):
        self.individuals.append(individual)
        self.individuals_by_id[individual.getId()] = individual
    def read_loci(self, locus_file):
        # print "read_loci()"
        fh = open(locus_file, "r")
        first = True
        count = 0L
        for line in fh.readlines():
            # Skip the header row
            if first:
                first = False
                continue
            tokens = self.whitespace.split(line.strip())
            assert len(tokens) == 3, "Locus file needs to have 3 columns."
            try:
                distance = float(tokens[2])
            except:
                distance = None
            if len(self.loci) > 0:
                prev_position = self.loci[-1].position
            else:
                prev_position = 0
            # print "Locus %s, distance %s, prev_position %s." \
            # % (tokens[0], distance, prev_position)
            locus = Locus(tokens[0], int(tokens[1]),
                    distance, prev_position, count)
            self.loci.append(locus)
            count += 1
        fh.close()
    def whitespace_split(self, line):
        """Accepts a list as a parameter and if it has only
        one element, splits it."""
        assert type(line) == list
        if len(line) == 1:
            return self.whitespace.split(line[0])
        else:
            return line
    def read_haplotypes(self, haplotype_file):
        # print "read_haplotypes()"
        fh = open(haplotype_file, "r")
        haplotypes = []
        # SNP IDs from the first line (header line)
        loci_snp_ids = self.whitespace_split([fh.readline(), ])[1:]
        count = 0L
        # Read the haplotypes
        while True:
            count += 1
            line1 = [fh.readline(),]
            if not line1[0]:
                break
            # Strip first column. The rest is haplotypes coded as
            # numbers
            line1 = self.whitespace_split(line1)
            id = line1[0]
            line1 = line1[1:]
            line2 = self.whitespace_split([fh.readline(),])[1:]
            indiv = Individual(id, count)
            try:
                assert(len(line1) == len(line2)), \
                    "Numbers of haplotypes need to match."
            except:
                print len(line1), len(line2)
                raise
            for i in range(len(line1)):
                genotypes = [line1[i], line2[i]]
                indiv.add_genotype(loci_snp_ids[i], genotypes)
            self.individuals.append(indiv)
            self.individuals_by_id[indiv.id] = indiv
            # For testing
            # if count > 5: break
        fh.close()

    def randomize_masking(self, percent_indivs, percent_loci, limit_loci):
        """Choose individuals and loci randomly."""
        print "Masking %s percent individuals and %s percent loci." \
                % (percent_indivs, percent_loci)
        assert len(self.individuals) > 0, "Individuals must be present"
        if not limit_loci or limit_loci > len(self.loci):
            limit_loci = len(self.loci)
        no_indivs_masked = len(self.individuals) * percent_indivs / 100
        if no_indivs_masked < 1:
            print "No individuals to mask."
            sys.exit(1)
        print "Masking %d individuals" % (no_indivs_masked)
        # Choose indices of masked individuals by shuffling a list
        # and trimming it.
        masked_indivs_list = range(len(self.individuals))
        random.shuffle(masked_indivs_list)
        masked_indivs_list = masked_indivs_list[:no_indivs_masked]
        for indiv_idx in masked_indivs_list:
            self.masked_indivs[self.individuals[indiv_idx].id] = 1

        # The same for loci.
        no_loci_masked = limit_loci * percent_loci / 100
        if no_loci_masked < 1:
            print "No loci to mask."
            sys.exit(1)
        mask_loci_idxs = range(limit_loci)
        random.shuffle(mask_loci_idxs)
        mask_loci_idxs = mask_loci_idxs[:no_loci_masked]
        for locus_idx in mask_loci_idxs:
            self.masked_loci[self.loci[locus_idx].snp_id] = self.loci[locus_idx]

        # Find out if there are any loci that have became monomorphic
        # because of masking. Those loci should be dropped.
        # Building a list of loci to write.
        self.loci_to_write = filter(
                lambda x: not x.is_monomorphic(self.get_unmasked_individuals()),
                self.loci[:limit_loci])
        count = 0L
        # print self.loci_to_write
        for locus in self.loci_to_write:
            locus.set_masked_idx(count)
            count += 1
        print "Number of loci to write: %s out of %s with a limit of %s" % \
                (len(self.loci_to_write), len(self.loci), limit_loci)

    def get_masked_loci_to_write(self):
        return filter(lambda x: self.masked_loci.has_key(x.snp_id), self.loci_to_write)

    def get_unmasked_loci(self):
        unmasked_loci = []
        for locus in self.loci:
            if not self.masked_loci.has_key(locus.snp_id):
                unmasked_loci.append(locus)
        return unmasked_loci

class HapmixmapFormatter:
    def format_haplotypes(self, indivs, loci, masked_loci, mask = False):
        ret = ""
        for indiv in indivs:
            for gn in range(2):
                line_list = ["%s_%s" % (indiv.id, gn),]
                line_list.extend(map(lambda x: indiv.get_genotype(x.snp_id)[gn], loci))
                ret += " ".join(line_list) + "\n"
        return ret
    def format_single_diplotype(self, g1, g2, snp_id, masked_loci, mask):
        assert type(mask) == bool
        assert type(snp_id) == str
        assert type(masked_loci) == dict
        if mask and masked_loci.has_key(snp_id):
            return '"0,0"'
        else:
            return '"%s,%s"' % (g1, g2)
    def format_diplotypes(self, indivs, loci, masked_loci, mask = False):
        ret = ""
        for indiv in indivs:
            line_list = ["%s" % (indiv.id),]
            line_list.extend(map(
                lambda x: self.format_single_diplotype(
                    indiv.get_genotype(x.snp_id)[0],
                    indiv.get_genotype(x.snp_id)[1],
                    x.snp_id,
                    masked_loci,
                    mask),
                loci))
            ret += "\t".join(line_list) + "\n"
        return ret
    def format_header(self, indivs, loci, masked_loci, indiv_id = "IndivId"):
        field_list = [indiv_id, ]
        field_list.extend(map(lambda x: x.snp_id, loci))
        return "\t".join(field_list) + "\n"

class FastPhaseFormatter:
    def write(self, file_name, indivs, loci, masked_loci):
        fh = open(file_name, "w")
        fh.write("%s\n" % len(indivs))
        fh.write("%s\n" % len(loci))
        fh.write("P")
        for locus in loci:
            fh.write(" %2.6f" % (locus.position * 1e-6))
        fh.write("\n")
        # fh.write("SSSSSSSSSSSSSSSSSSS\n") # Funny format, isn't it?
        for indiv in indivs:
            fh.write("# %s\n" % indiv.id)
            # Two lines with data per individual
            for gn in range(2):
                for locus in loci:
                    if masked_loci.has_key(locus.snp_id):
                        fh.write("?")
                    else:
                        fh.write(indiv.get_genotype(locus.snp_id)[gn])
                fh.write("\n")
        fh.close()

class HapmixmapIndexFormatter:
    def write(self, file_name, indivs, loci, masked_loci, unmasked_indivs):
        assert type(masked_loci) == dict
        for locus_id in masked_loci.keys():
            assert masked_loci[locus_id].__class__.__name__ == 'Locus'
            if masked_loci[locus_id] in loci:
                assert not masked_loci[locus_id].masked_idx is None
        fh = open(file_name, "w")
        fh.write("maskedindivs = ")
        # The following code would be used if the masked and unmasked
        # were in their original position.
        #
        # Extract the indices from individuals
        # indivs_idxs = map(lambda x: x.idx, indivs)
        # Sort the indices numerically
        # indivs_idxs.sort()
        #
        # What the program is going to use now, is that the masked
        # individuals are the last ones; in the set. The indices refer
        # to 1-based positions in the merged (haploid + diploid) file.
        # Number of gametes equals number of individuals * 2
        masked_idx_start = len(unmasked_indivs) * 2
        indivs_idxs = range(masked_idx_start, masked_idx_start + len(indivs))
        
        # Convert them to strings
        # WARNING: hapmixmap expects 1-based indices (hence +1)
        indivs_idxs = map(lambda x: str(x + 1), indivs_idxs)
        fh.write(" ".join(indivs_idxs))
        fh.write("\n")
        fh.write("maskedloci = ")
        # Some masked loci may be have been excluded (have become
        # monomorphic) and don't have masked_idx. These are filtered
        # out.
        masked_loci_keys = filter(lambda x: masked_loci[x] in loci, masked_loci.keys())
        loci_idxs = map(lambda x: masked_loci[x].masked_idx, masked_loci_keys)
        loci_idxs.sort()
        # WARNING: hapmixmap expects 1-based indices (hence +1)
        fh.write(" ".join(map(lambda x: str(x + 1), loci_idxs)))
        fh.write("\n")
        fh.close()

# R-project format, can be read by `dget' command.
#
# structure(list(rsXXXX = c("2,1", ...), ...),
# .Names = c("rsXXXXX", ...), row.names = c("NA0001", ...),
# class = "data.frame")
#
class DgetFormatter:
    def format(self, indivs, loci, masked_loci):
        assert type(loci) == list
        assert loci[0].__class__.__name__ == 'Locus'
        # TODO: continue from here, R formatter
        ret = ""
        ret += "structure(list(\n"
        # ret += "%s = c(" % (locus.snp_id)

        # Build diploid data structure
        loci_list = []
        for locus in loci:
            locus_dd = {}
            locus_dd['genotypes'] = []
            locus_dd['snp_id'] = locus.snp_id
            for indiv in indivs:
                locus_dd['genotypes'].append('"%s,%s"' % (
                        indiv.get_genotype(locus.snp_id)[0],
                        indiv.get_genotype(locus.snp_id)[1],))
            loci_list.append(locus_dd)
        
        # Format the data (R-specific)
        ret += ",\n".join(map(lambda x:
                ("%s = c(" % (x['snp_id']))
                + ", ".join(x['genotypes'])
                + ")",
                loci_list))
        ret += "),\n"

        # .Names = c("rsXXXXX", ...), row.names = c("NA0001", ...),
        ret += ".Names = c("
        ret += ", ".join(map(lambda x: '"%s"' % x.snp_id, loci))
        ret += "),\n"
        ret += "row.names = c("
        ret += ", ".join(map(lambda x: '"%s"' % x.id, indivs))
        ret += "),\n"
        ret += 'class = "data.frame")'

        return ret

class LocusFileFormatter:
    # When writing a lost of loci, the distances need to be recalculated
    def format(self, loci):
        ret = '"SNPid" "NumAlleles"    "DistanceinMb"\n'
        prev_pos = 0
        count = 0L
        for locus in loci:
            if 0 == count:
                distance = "#"
            else:
                distance = "%0.8f" % (1e-6 * (locus.position - prev_pos), )
            locus_list = [locus.snp_id, "2", distance]
            ret += "\t".join(locus_list) + "\n"
            prev_pos = locus.position
            count += 1
        return ret

class DataMaskingMediator:
    """Currently only handles paths."""
    def __init__(self, out_dir):
        self.output_dir = out_dir
    def get_output_dir(self):
        return self.output_dir
    def get_output_filename(self, fn):
        return os.path.join(self.output_dir, fn)


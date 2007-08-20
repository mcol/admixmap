#!/usr/bin/python
# -*- coding: UTF-8 -*-
#
# Classes for manipulating genetic data.
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
import time
import gc
import itertools
from threading import Thread

class ProgressMeter(object):
    """Monitor progress of something."""
    def __init__(self, fh, total, every = 1):
        self.every = every
        self.counter = 0L
        self.fh = fh
        self.total = total
    def step(self):
        self.counter += 1
        if (self.counter % self.every) == 0:
            self.fh.write("\r%d / %d (%02.2f%%)" % (self.counter, self.total, 100.0 * self.counter / self.total))
            self.fh.flush()
    def finalize(self):
        self.fh.write("\n")
        self.fh.flush()

class TimeReporter(object):
    def __init__(self):
        self.starttime = time.clock()
        self.last = time.clock()
    def checkpoint(self, msg):
        newtime = time.clock()
        print "%s" % msg, newtime - self.last
        self.last = newtime

class Locus(object):
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
    def getSnpId(self):
        return self.snp_id
    def get_coding(self, conn):
        curs = conn.cursor()
        qry = """
        SELECT
            l.snp_id,
            l.chromosome,
            l.position,
            l.allele1,
            l.allele2
        FROM
            locus AS l
        WHERE
            l.snp_id = %s
        ;
        """
        curs.execute(qry, (self.getSnpId(),));
        locus = curs.fetchone()
        return [locus[3], locus[4]]
    def __repr__(self):
        return "<Locus: %s>" % self.snp_id

class Individual(object):
    genotype_cache = {}
    snp_id_cache = {}
    def __init__(self, id, idx):
        self.idx = idx
        self.genotypes = {}
        self.id = id
        # self.genotypes_list = []
        # self.genotypes_ready = False
        # print "Individual(%s)" % repr(id)
    def add_genotype(self, locus_id, genotype):
        # print "add_genotype(%s, %s)" % (locus_id, genotype)
        # self.genotypes_list.append((locus_id, genotype))
        genotype2 = self.genotype_cache.setdefault(genotype, genotype)
        locus_id2 = self.snp_id_cache.setdefault(locus_id, locus_id)
        self.genotypes[locus_id2] = genotype2
    # def prepare_genotypes(self):
        # if not self.genotypes_ready:
        #     for locus_id, genotype in self.genotypes_list:
        #         self.genotypes[locus_id] = genotype
        #     self.genotypes_ready = True
    def getId(self):
        return self.id
    def get_genotype(self, snp_id):
        # self.prepare_genotypes()
        if self.genotypes.has_key(snp_id):
            return self.genotypes[snp_id]
        else:
            # Perhaps throw an exception?
            print "Snp %s not found in individual %s." % \
                    (snp_id, self.id)
            sys.exit(1)

class UserData(object):
    def __init__(self):
        self.individuals = []
        self.individuals_by_id = {}
        self.loci = []
        self.loci_by_id = {}
        self.whitespace = re.compile("\s+")
        self.masked_loci = {}
        self.masked_indivs = {}
        self.loci_to_write = []
        self.loci_count = 0
        self.tr = TimeReporter()
    def remove_individuals(self):
        self.individuals = []
        self.individuals_by_id = {}
    def get_unmasked_individuals(self):
        return filter(lambda x: not self.masked_indivs.has_key(x.id),
                self.individuals)
    def get_masked_individuals(self):
        return filter(lambda x: self.masked_indivs.has_key(x.id),
                self.individuals)
    def add_individual(self, individual):
        self.individuals.append(individual)
        self.individuals_by_id[individual.getId()] = individual
    def read_hapmap_legend(self, legend_file):
        fh = open(legend_file, "r")
        prev_pos = 0
        lines = fh.readlines()
        # print "Reading", len(lines), "lines of locus file, containing",
        # print len(lines) - 1, "loci."
        lines = lines[1:]
        for line in lines:
            tokens = self.whitespace.split(line.strip())
            try:
                assert len(tokens) == 4, "Four tokens needed"
            except:
                print tokens
                raise

            locus = Locus(tokens[0], 2, 1e-6 * (prev_pos - int(tokens[1])),
                    1e-6 * prev_pos, self.loci_count)
            self.addLocus(locus)
            self.loci_count += 1
    def read_loci(self, locus_file):
        """Read loci in hapmixmap format."""
        # print "read_loci()"
        fh = open(locus_file, "r")
        first = True
        # count = 0L
        for line in fh.readlines():
            # Skip the header row
            if first:
                first = False
                continue
            tokens = self.whitespace.split(line.strip())
            try:
                assert len(tokens) == 3, "Locus file needs to have 3 columns."
            except:
                print tokens
                raise
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
                    distance, prev_position, self.loci_count)
            self.addLocus(locus)
            self.loci_count += 1
        fh.close()
    def addLocus(self, l):
        self.loci.append(l)
    def whitespace_split(self, line):
        """Accepts a list as a parameter and if it has only
        one element, splits it."""
        assert type(line) == list
        if len(line) == 1:
            return self.whitespace.split(line[0].strip())
        else:
            return line
    def read_haplotypes(self, haplotype_file):
        """Read haplotypes in hapmixmap format."""
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
            self.addIndividualByLines(id, count, line1, line2,
                    loci_snp_ids)
            # For testing
            # if count > 5: break
        fh.close()
    def read_hapmap_haplotypes(self, haplotypes_fn, sample_fn, batch = []):
        """Read haplotypes in hapmap format."""
        # Hapmap encodes genotypes as 0 and 1 while Hapmixmap encodes
        # them as 1 and 2.
        gt_coding = {
                '0': '1',
                '1': '2', }
        # gc.disable()
        indiv_time = TimeReporter()
        sample_fh = open(sample_fn, "r")
        sample = map(lambda x: self.whitespace_split([x]), sample_fh.readlines())
        assert len(sample) > 0
        fh = open(haplotypes_fn, "r")
        loci_snp_ids = [x.getSnpId() for x in self.loci]
        assert len(loci_snp_ids) > 0, "SNP list empty"
        # print "SNP id list contains", len(loci_snp_ids), "values."
        haplotypes = []
        count = 0L
        while True:
            # self.tr.checkpoint("indiv-start: ")
            line1 = [fh.readline(),]
            if not line1[0]: break
            raw_line2 = fh.readline()
            if count in batch or not batch:
                line1 = self.whitespace_split(line1)
                line2 = self.whitespace_split([raw_line2,])
                line1 = [gt_coding[x] for x in line1]
                line2 = [gt_coding[x] for x in line2]
                self.addIndividualByLines(sample[count][0], count,
                        line1, line2, loci_snp_ids)
                # self.tr.checkpoint("individual-added: ")
                if not (count % 10):
                    print "Individual no:", count,
                    indiv_time.checkpoint("time: ")
            # self.tr.checkpoint("lines splitted: ")
            # print "Looking for individual id. Count is", count,
            # print "and the sample[count] is", sample[count]
            count += 1
            # For testing
            # if count > 5: break
            # if 0 == count % 10:
            #     gc.collect()
        # print "Enabling the garbage collector."
        # gc.enable()

    def addIndividualByLines(self, id, no, line1, line2, loci_snp_ids):
        indiv = Individual(id, no)
        try:
            assert len(line1) == len(line2), \
                "Numbers of haplotypes need to match."
        except:
            print len(line1), len(line2)
            raise
        try:
            assert len(line1) == len(loci_snp_ids), \
                    "Number of SNP ids needs to match line elements"
        except:
            print "line1:", len(line1), "loci_snp_ids:", len(loci_snp_ids)
            print line1[:5], loci_snp_ids[:5]
            print line1
            raise
        # self.tr.checkpoint("starting to add genotypes: ")
        # print "Adding", len(line1), "genotypes."
        for i in range(len(line1)):
            indiv.add_genotype(loci_snp_ids[i], (line1[i], line2[i]))
        # self.tr.checkpoint("genotypes added, appending individual: ")
        self.individuals.append(indiv)
        self.individuals_by_id[indiv.id] = indiv
        del indiv
        # self.tr.checkpoint("individual appended: ")

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

class HapmixmapFormatter(object):
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

class FastPhaseFormatter(object):
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

class HapmixmapIndexFormatter(object):
    def write(self, file_name, indivs, loci, masked_loci, unmasked_indivs):
        assert type(masked_loci) == dict
        for locus_id in masked_loci.keys():
            assert isistance(masked_loci[locus_id], Locus), "Locus needed"
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
class DgetFormatter(object):
    def format(self, indivs, loci, masked_loci):
        assert type(loci) == list
        assert loci[0].__class__.__name__ == 'Locus'
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

class SqlFormatter(object):
    def __init__(self, conn, loci_numbers):
        import psycopg
        self.conn = conn
        self.curs = self.conn.cursor()
        self.loci_numbers = loci_numbers
    def write_header(self, fh):
        fh.write("COPY observed_genotype (indiv_id, locus_id, allele1, allele2) FROM stdin;\n")
    def write_footer(self, fh):
        fh.write("\\.\n")
    def write(self, fh, indivs, loci, data_set):
        assert type(loci) == list, "Loci must be a list"
        assert len(loci) > 0, "Loci list must be non-empty"
        assert isinstance(loci[0], Locus), "Locus class elements needed"
        # Check if the data set exists
        loci_numbers = self.loci_numbers
        self.curs.execute("SELECT * FROM data_set WHERE data_set_name = %s;",
                (data_set, ))
        if self.curs.rowcount == 0:
            self.curs.execute("INSERT INTO data_set (data_set_name) VALUES (%s);",
                    (data_set, ))
            self.conn.commit()
        for indiv in indivs:
            self.curs.execute("SELECT id FROM individual WHERE char_id = %s AND data_set_name = %s;", (indiv.getId(), data_set))
            only_row = self.curs.fetchone()
            if self.curs.rowcount == 0:
                print "inserting individual", indiv.getId()
                self.curs.execute("""
                INSERT INTO individual (char_id, data_set_name)
                VALUES (%s, %s); """, ( indiv.getId(), data_set))
                self.conn.commit()
                # get the new individual's ID
                self.curs.execute("SELECT id FROM individual WHERE char_id = %s;", (indiv.getId(), ))
                only_row = self.curs.fetchone()
            else:
                print "processing data for", indiv.getId()
            indiv_no = only_row[0]
            assert(indiv_no)
            # self.curs.free()
            # progressMeter = ProgressMeter(sys.stdout, len(loci), 1000)
            for locus in loci:
                # sys.stdout.write("locus %s\n" % locus.getSnpId())
                # sys.stdout.flush()
                # progressMeter.step()
                snp_id = locus.getSnpId()
                if not loci_numbers.has_key(snp_id):
                    self.curs.execute("SELECT id FROM locus WHERE snp_id = %s;", (locus.getSnpId(), ))
                    loci_numbers[snp_id] = self.curs.fetchone()[0]
                    # self.curs.free()
                fh.write("""%(indiv_no)s\t%(locus_id)s\t%(allele1)s\t%(allele2)s \n""" % { 'indiv_no': indiv_no,
        'locus_id': loci_numbers[snp_id],
        'allele1': indiv.get_genotype(snp_id)[0],
        'allele2': indiv.get_genotype(snp_id)[1], })
            # progressMeter.finalize()

class LocusFileFormatter(object):
    # When writing a lost of loci, the distances need to be recalculated
    def format(self, loci):
        ret = '"SNPid" "NumAlleles"    "DistanceinMb"\n'
        prev_pos = 0
        count = 0L
        for locus in loci:
            if 0 == count:
                distance = "#"
            else:
                distance = "%01.8f" % (1e-6 * (prev_pos - locus.position), )
            locus_list = [locus.snp_id, "2", distance]
            ret += "\t".join(locus_list) + "\n"
            prev_pos = locus.position
            count += 1
        return ret

class DataMaskingMediator(object):
    """Currently only handles paths."""
    def __init__(self, out_dir):
        self.output_dir = out_dir
    def get_output_dir(self):
        return self.output_dir
    def get_output_filename(self, fn):
        return os.path.join(self.output_dir, fn)

class SqlConverter(object):
    def __init__(self, conn, haplo_fn, locus_fn, sam_fn, loci_numbers):
        self.conn = conn
        self.haplotypes_fn = haplo_fn
        self.locus_fn = locus_fn
        self.sample_fn = sam_fn
        self.ud = UserData()
        self.curs = self.conn.cursor()
        self.loci_numbers = loci_numbers

    def write(self, fh, data_set_name, max_indivs = 21):
        tot_indivs = len(open(self.haplotypes_fn, "r").readlines()) / 2
        batches = make_batches(tot_indivs, max_indivs)
        print "Reading hapmap loci"
        self.ud.read_hapmap_legend(self.locus_fn)
        sdf = SqlFormatter(self.conn, self.loci_numbers)
        sdf.write_header(fh)
        for batch in batches:
            print "Reading hapmap haplotypes for", batch
            self.ud.remove_individuals()
            self.ud.read_hapmap_haplotypes(
                    self.haplotypes_fn, self.sample_fn, batch)
            print "Calling SqlFormatter::write()", batch
            sdf.write(fh, self.ud.individuals, self.ud.loci, data_set_name)
            gc.collect()
        sdf.write_footer(fh)
        del sdf

class SqlInserter(SqlConverter):
    def write(self, data_set_name):
        print "Reading hapmap loci"
        self.ud.read_hapmap_legend(self.locus_fn)
        self.ud.read_hapmap_haplotypes(
                self.haplotypes_fn, self.sample_fn) # , range(8))
        loci_numbers = self.loci_numbers
        pm = ProgressMeter(sys.stdout, len(self.ud.loci), every = 500)
        count = 0L
        indiv_ids = ""
        for indiv in self.ud.individuals:
            indiv_ids += "%s\n" % indiv.getId()
        qry = "UPDATE data_set SET indiv_ids = %s WHERE data_set_name = %s"
        self.curs.execute(qry, [indiv_ids, data_set_name])
        for locus in self.ud.loci:
            gt_set = set([])
            mm = True
            snp_id = locus.getSnpId()

            # DEBUG
            # if snp_id == 'rs8139476':
            #     print "got rs8139476 for %s" % (data_set_name, )

            # Create a text object with all the genotype values for
            # this individual.
            txt = ""
            for indiv in self.ud.individuals:
                gt = indiv.get_genotype(snp_id)
                txt += "%s\n%s\n" % (gt[0], gt[1])
                gt_set.add(gt[0])
                gt_set.add(gt[1])
            # print "locus: %s, txt: %s" % (snp_id, txt)
            if len(gt_set) > 1:
                mm = False
            qry = """
            SELECT lp.id FROM locus_population AS lp
            INNER JOIN locus AS l ON (lp.locus_id = l.id)
            WHERE l.snp_id = %s AND lp.data_set_name = %s;
            """
            self.curs.execute(qry, (snp_id, data_set_name))
            row = self.curs.fetchone()

            # if snp_id == 'rs8139476':
            #     if not row:
            #         print "row NOT found for rs8139476, %s" % (data_set_name, )
            #         print qry, (snp_id, data_set_name)
            #     else :
            #         print "row FOUND for rs8139476, %s" % (data_set_name, )

            if not row:
                qry = """
    INSERT INTO locus_population (data_set_name, locus_id, monomorphic, indivs)
    VALUES (%s, (SELECT id from locus where snp_id = %s), %s, %s);
                """

                # if snp_id == 'rs8139476':
                #     print qry, (data_set_name, snp_id, mm, 'txt')

                self.curs.execute(qry, (data_set_name, snp_id, mm, txt))
            pm.step()
            if not (count % 500):
                self.conn.commit()
            count += 1
        self.conn.commit()
        pm.finalize()

def make_batches(no_elems, per_batch):
    """Returns a list of lists with indexes."""
    no_batches = no_elems / per_batch
    rest = no_elems % per_batch
    batches = []
    for batch in range(no_batches):
        batches.append([x + batch * per_batch for x in range(per_batch)])
    batches.append([x + no_batches * per_batch for x in range(rest)])
    return batches

class Hapmap(object):
    def __init__(self):
        pass
    def get_legend_fn(self, chromosome, pop):
        return "genotypes-%s-%s-r21-nr-fwd-legend.txt" \
                % (chromosome, pop)
    def get_haplotypes_fn(self, chromosome, pop):
        return "genotypes-%s-%s-r21-nr-fwd-phased" \
                % (chromosome, pop)
    def get_sample_fn(self, chromosome, pop):
        return "genotypes-%s-%s-r21-nr-fwd-sample.txt" \
                % (chromosome, pop)
    def get_populations(self):
        return ['ceu', 'jpt-chb', 'yri']
    def get_chromosomes(self):
        return ["chr%s" % x for x in range(1, 23)]

class ColumnFilesConverter(object):
    def __init__(self, haplo_fn, locus_fn, sam_fn, target_dir):
        self.haplotypes_fn = haplo_fn
        self.locus_fn = locus_fn
        self.sample_fn = sam_fn
        self.ud = UserData()
        self.target_dir = target_dir
    def write(self, data_set_name, max_indivs = 21):
        tot_indivs = len(open(self.haplotypes_fn, "r").readlines()) / 2
        batches = make_batches(tot_indivs, max_indivs)
        print "Reading hapmap loci"
        self.ud.read_hapmap_legend(self.locus_fn)
        for batch in batches:
            print "Reading hapmap haplotypes for", batch
            self.ud.remove_individuals()
            self.ud.read_hapmap_haplotypes(
                    self.haplotypes_fn, self.sample_fn, batch)
            self.writeColumns(data_set_name)
            gc.collect()
    def writeColumns(self, data_set_name):
        print "ColumnFilesConverter::writeColumns(%s)" % data_set_name
        total_loci = len(self.ud.loci)
        count = 0L
        for locus in self.ud.loci:
            snp_id = locus.getSnpId()
            column_fn = self.getColumnFileName(data_set_name, snp_id)
            # print column_fn
            # Check if 1-st level dir exists
            cfpl = self.getColumnFilePathList(snp_id, data_set_name)
            assert len(cfpl) == 5, "Path needs to consists of 4 elements"
            # Check all the dirs, if they not exist, create them.
            for i in range(2, 5):
                check_dir = os.path.join(*cfpl[:i])
                if not os.path.isdir(check_dir):
                    os.mkdir(check_dir)
            fh = open(os.path.join(*cfpl), "a")
            for indiv in self.ud.individuals:
                fh.write("%s\n" % indiv.get_genotype(snp_id)[0])
                fh.write("%s\n" % indiv.get_genotype(snp_id)[1])
            fh.close()
            count += 1
            if count % 10000 == 0:
                print "loci count:", count, "of", total_loci,
                print " %02.2f%%" % (100.0 * float(count) / float(total_loci))
    def getColumnFileName(self, data_set_name, snp_id):
        return os.path.join(*self.getColumnFilePathList(snp_id, data_set_name))
    def getColumnFilePathList(self, snp_id, data_set_name):
        path_list = [self.target_dir, data_set_name]
        path_list.extend(self.getSnpFileBase(snp_id))
        return path_list
    def getSnpFileBase(self, snp_id):
        return [self.getFirstLevelDir(snp_id),
            self.getSecondLevelDir(snp_id),
            self.getBaseColumnFile(snp_id)]
    def getFirstLevelDir(self, snp_id):
        return snp_id[:3]
    def getSecondLevelDir(self, snp_id):
        return snp_id[:4]
    def getBaseColumnFile(self, snp_id):
        return "%s.txt" % snp_id

class ImputeFormatter(object):
    gt_map = {
            '1': '0',
            '2': '1', }
    dip_probs = {
            '"0,0"': '0 0 0',
            '"1,1"': '1 0 0',
            '"1,2"': '0 1 0',
            '"2,1"': '0 1 0',
            '"2,2"': '0 0 1', }
    def write_haplo(self, fh, indivs, loci):
        """
        Rows: SNPs
        Columns: Individuals (gametes?)

        0 0 0 0 1 1 0 0 0 1 0 ...
        1 1 1 0 0 0 1 1 1 1 1 ...
        """
        assert isinstance(fh, file)
        assert isinstance(indivs, list)
        assert len(indivs) > 0, "List is empty"
        assert isinstance(indivs[0], Individual)
        for locus in loci:
            snp_id = locus.getSnpId()
            line = " ".join([
                "%s %s" % (
                self.gt_map[i.get_genotype(snp_id)[0]],
                self.gt_map[i.get_genotype(snp_id)[1]])
                for i in indivs])
            fh.write(line)
            fh.write("\n")
    def write_geno(self, fh, indivs, loci, conn):
        """
SNP_A-2173923 rs9608606 26076638 A G 0 0 1 1 0 0 0 1 0
SNP_A-4294550 rs2179099 26079328 C T 0 0 1 0 0 1 0 1 0
        """
        tmpl = "%(snp_id)s %(snp_id)s %(position)s %(g1)s %(g2)s "
        for locus in loci:
            snp_id = locus.getSnpId()
            coding = locus.get_coding(conn)
            line = tmpl % {
                'snp_id': snp_id,
                'position': locus.position,
                'g1': coding[0],
                'g2': coding[1], }
            line += " ".join([
                self.dip_probs[i.get_genotype(snp_id)]
                for i in indivs])
            fh.write(line)
            fh.write("\n")
    def write_legend(self, fh, loci, conn):
        """
rs position X0 X1
rs9613418 26032219 A G
rs13058071 26032867 A C
rs528686 26033406 A G
        """
        fh.write("rs position X0 X1\n")
        for locus in loci:
            coding = locus.get_coding(conn)
            fh.write("%s %s %s %s\n" % (locus.snp_id, locus.position,
                coding[0], coding[1]))


class DiploidUserData(UserData):
    """Handles diploid data."""
    def read_diplotypes(self, fh):
        I_CUTE = True
        assert isinstance(fh, file)
        loci_snp_ids_file = self.whitespace_split([fh.readline(), ])[1:]
        loci_snp_ids = [x.snp_id for x in self.loci]
        assert loci_snp_ids_file == loci_snp_ids
        assert len(loci_snp_ids) > 0, "SNP list empty"
        count = 0L
        while I_CUTE: # LOL
            line1 = [fh.readline(),]
            if not line1[0]: break
            line1 = self.whitespace_split(line1)
            # First column is individual id
            id = line1[0]
            line1 = line1[1:]
            self.addIndividualByLines("IND_%s" % count, count,
                    line1, loci_snp_ids)
            count += 1
    def addIndividualByLines(self, id, no, line1, loci_snp_ids):
        indiv = DiploidIndividual(id, no)
        for i in range(len(line1)):
            try:
                indiv.add_genotype(loci_snp_ids[i], line1[i])
            except IndexError, e:
                print "i =", i,
                print "len(loci_snp_ids) =", len(loci_snp_ids),
                print "len(line1) =", len(line1),
                print
                print e
                raise
        self.individuals.append(indiv)
        self.individuals_by_id[indiv.id] = indiv

class DiploidIndividual(Individual):
    """Individual with diploid data."""
    pass

def group_list(lst, n):
    """Changes [1, 2, 3, 4, ...] into [[1, 2], [3, 4], ...]"""
    return map(lambda x: tuple(lst[x:x + n]),
            itertools.islice(range(len(lst)), 0, len(lst), n))

class RArray(object):
    columns = 3
    def __init__(self, data, dims, labels):
        """
        data is a list
        dims is a list of dimension lengths
        labels is a list of lists of strings
        """

        # It can be an iterator
        # assert isinstance(data, list)
        assert isinstance(dims, list)
        assert isinstance(labels, list)
        assert len(dims) == len(labels), "dims: %s, labels: %s" \
                % (len(dims), len(labels))
        self.data = data
        self.dims = dims
        self.labels = labels
    def __cols(self, cols = None):
        if not cols:
            return self.columns
        else:
            return cols
    def __format_columns(self, lst, cols = None, quote = False):
        cols = self.__cols(cols)
        if quote:
            quote_s = '"'
        else:
            quote_s = ""
        lines = [lst[i:i + cols] for i in xrange(0, len(lst), cols)]
        lines = [", ".join(["%s%s%s" % (quote_s, y, quote_s) for y in x]) for x in lines]
        return ",\n".join(lines)
    def __format_columns_lol(self, lst, cols = None, quote = False):
        cols = self.__cols(cols)
        return ",\n".join(["c(%s)" % y for y in [self.__format_columns(x, quote = quote) for x in lst]])
    def format_data(self, d):
        # print "RArray::format_data()"
        return self.__format_columns(d)

    def format(self):
        r = "# vim:set ft=r:\n"
        r += "structure(.Data = c(\n"
        # Formatting values, using specified number of self.columns
        # [[0, 1, 2], [3]]
        # ['0, 1, 2', '3']
        r += self.format_data(self.data)
        r += "),\n"
        # Dimensions
        r += "\t.Dim = c(" + self.__format_columns(self.dims, quote = True) + "),\n" 
        # Dimension labels
        r += "\t.Dimnames = list(\n" + self.__format_columns_lol(self.labels, quote = True) + "))\n" 
        return r

class LociPriors(object):
    def __init__(self, loci, pop, conn):
        self.conn = conn
        self.loci = loci
        self.pop = pop
        self.snp_ids = [x.snp_id for x in self.loci]
        self.by_snp_id = {}
        for locus in self.loci:
            self.by_snp_id[locus.snp_id] = {'locus': locus, }
    def format_dget(self):
        curs = self.conn.cursor()
        for snp_id in self.snp_ids:
            qry = """
SELECT
    snp_id, lp.id, gametes, a1, a2
FROM
    locus AS l
    INNER JOIN locus_population AS lp ON (l.id = lp.locus_id)
WHERE
    lp.data_set_name = %s
    AND
    snp_id = %s
;
            """
            curs.execute(qry, [self.pop, snp_id])
            row = curs.fetchone()
            assert row, "No data in the database, pop: %s, snp_id: %s" % (self.pop, snp_id)
            assert row[2] == row[3] + row[4], "Wrong values: %s != %s + %s" % row[2:5]
            p = float(row[3]) / row[2]
            q = float(row[4]) / row[2]
            tmp_priors = (p * p, 2 * p * q, q * q)
            assert sum(tmp_priors) - 1.0 < 1e-8, \
                    "Priors need to sum to 1.0: sum(%s) = %s" % (tmp_priors, sum(tmp_priors))
            self.by_snp_id[snp_id]['priors'] = tmp_priors
        # for snp_id in self.snp_ids:
        #     print self.by_snp_id[snp_id]
        data = []
        data.extend([self.by_snp_id[x]['priors'][0] for x in self.snp_ids])
        data.extend([self.by_snp_id[x]['priors'][1] for x in self.snp_ids])
        data.extend([self.by_snp_id[x]['priors'][2] for x in self.snp_ids])
        dims = [len(self.snp_ids), 3]
        labels = []
        # row labels
        labels.append(self.snp_ids)
        # column labels (0, 1, 2)
        labels.append(range(3))
        array = RArray(data, dims, labels)
        return array.format()

class HapmixmapMaskingIndex(object):
    def __init__(self, fn):
        self.fn = fn
        self.idx = set()
        self.idx_list = []
        for line in open(fn, "r").readlines():
            if re.search(r'maskedloci', line):
                self.process_line(line)
    def process_line(self, line):
        numbers = re.findall(r'[0-9]+', line.split("=")[1])
        for num in numbers:
            n = int(num)
            self.idx.add(n)
            self.idx_list.append(n)
    def filter_loci(self, loci):
		"""Takes a list of loci and returns a subset of that list."""
		# self.idx is 1-based
		return map(lambda x: loci[x - 1], self.idx)

if __name__ == '__main__':
    print "This is a library file."
    sys.exit(1)

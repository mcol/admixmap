from gzip import open as gzip_open
import re
import os.path
import sys

class ProgressIndicator:
    def __init__(self, target_no):
        self.counter = 0L
        self.target_no = target_no
        self.granularity = int(target_no / 1000)
    def update(self, msg = ""):
        self.counter += 1L
        if (self.counter % self.granularity) == 0:
            sys.stdout.write("\r%02.1f%% %s" % ((float(self.counter) / float(self.target_no)), msg))

class HapmapChromosomeLegend:
    def __init__(self, full_path):
        file_name = os.path.split(full_path)[-1]
        print "HapmapChromosomeLegend(%s)" % file_name
        if (not re.search(r'^genotypes.*\.txt\.gz$', file_name)):
            print "Can't initialize."
            return
        # genotypes_chr15_YRI_r21_nr_fwd_legend.txt.gz
        # legend_name_re = re.compile(r'genotypes_(?P<chr>)_(?P<pop>)_.*\.txt\.gz')
        legend_name_re = re.compile(r'genotypes_chr(?P<chr>[^_]+)_(?P<pop>[^_]+)_(?P<r>[^_]+)_nr_fwd(_par2)?_legend\.txt\.gz')
        srch = re.search(legend_name_re, file_name)
        if not srch:
            print "Can't parse the file name."
            return
        if srch:
            print srch.group('chr'), srch.group('pop'), srch.group('r')
        # self.file = gzip_open(file_name, "r")
        self.file = gzip_open(full_path, "r")
        self.rs = {}
    def read_rs_dict(self):
        self.file.seek(0)
        # header row
        no_lines = len(self.file.readlines()) - 1
        p_ind = ProgressIndicator(no_lines)
        self.file.seek(0)
        # Strip the header line.
        self.file.readline()
        while True:
            line_list = self.file.readline().strip().split('\t')
            if not line_list: break
            # self.rs[line_list[0]] = line_list[1:]
            # print "'%s': %s," % (line_list[0], line_list[1:])
            p_ind.update()
    def test(self):
        self.read_rs_dict()
        for key in self.rs.keys()[:10]:
            print key, self.rs[key]


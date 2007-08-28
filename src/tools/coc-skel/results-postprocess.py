#!/usr/bin/python
# -*- coding: UTF-8 -*-
# 
# HTML generator for Hapmixmap results
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
import os, os.path, shutil, stat
import re
from optparse import OptionParser
from hapmixmap import TextTable, HapmixmapLog, PostScriptImage, TextFile
from hapmixmap import slugify

# For compatibility with Python 2.3
if sys.version_info[0] <= 2 and sys.version_info[1] <= 3:
    from sets import Set as set

# Unavailable in Python 2.3
# from subprocess import Popen, PIPE

# p = Popen(file, shell=True,stdout=PIPE,stderr=PIPE)
# out = string.join(p.stdout.readlines())
# outerr = string.join(p.stderr.readlines())

class ResultsDirContainer(object):
    def __init__(self, out_dir):
        self.results = set()
        self.out_dir = out_dir
        self.ordered = dict()
    def add_dir(self, d):
        new_dir = ResultsDir(d, self.out_dir)
        self.results.add(new_dir)
        self.ordered[new_dir.get_order_key()] = new_dir
    def generate_html(self):
        for rdir in self.results:
            rdir.preprocess()
        html = self.get_html()
        open(os.path.join(self.out_dir, "index.html"), "w").write(html)
        shutil.copyfile("template/hapmixmap.css", os.path.join(self.out_dir, "hapmixmap.css"))
    def get_html_toc(self):
        r = "<ul id=\"toc\">\n"
        r += "<li><a href=\"#trace-plots\">Trace plots</a></li>\n"
        for rdir in self.get_ordered_results():
            r += "<li><a href=\"%s\">%s</a></li>\n" % \
                    ("#%s" % slugify(rdir.path), rdir.get_pretty_name(), )
        r += "</ul>\n"
        return r
    def get_ordered_results(self):
        keys = self.ordered.keys()
        keys.sort()
        return [self.ordered[key] for key in keys]
    def get_html(self):
        idx_content = ""
        tmpl = open("template/index.html").read()
        for rdir in self.get_ordered_results():
            idx_content += rdir.get_html(self.out_dir)
        # Gather all the latest energy trace plots and display them
        # together.
        idx_content += self.get_energy_plots()
        return tmpl % {
            'toc': self.get_html_toc(),
            'index_table': idx_content,
            'date': 'date TODO', }
    def get_energy_plots(self):
        r = "<a name=\"trace-plots\"></a><h2>Energy trace plots</h2>\n"
        r += "<ul class=\"results_charts\">\n"
        for rdir in self.get_ordered_results():
            if rdir.images_by_name.has_key('EnergyTracePlot'):
                r += """
                <li><div style="float: left;"><a href="#%s" title="%s details">%s</a> <br />%s</div></li>
                """ % (rdir.get_slug(), rdir.get_pretty_name(), rdir.get_pretty_name() , rdir.images_by_name['EnergyTracePlot'].get_html())
        r += "</ul>\n"
        r += "<div style=\"clear: both;\"></div>\n"
        return r

class ResultsDir(object):
    txt_tables = ('PosteriorQuantiles.txt', )
    txt_file_names = ('logfile.txt', 'args.txt', )
    def __init__(self, d, out_dir):
        # print "ResultsDir(%s)" % repr(d)
        self.path = d
        self.images = []
        self.out_dir = out_dir
        self.log = HapmixmapLog(os.path.join(self.path, "logfile.txt"))
        self.txt_files = []
        for fn in self.txt_file_names:
            self.txt_files.append(TextFile(self.path, fn, self.out_dir))
        self.images_by_name = {}
    def preprocess(self):
        # self.run_r()
        self.convert_images()
        self.copy_files()
    def run_r(self):
        # Run R script
        command = "RESULTSDIR=\"%s\" R CMD BATCH --no-save --no-restore test/admixmap/AdmixmapOutput.R" % self.path
        print command; out = os.popen(command).readlines()
    def convert_images(self):
        # preprocess graphics, 
        for fn in os.listdir(self.path):
            base_fn = fn[:-3]
            if fn[-3:] == ".ps":
                complete_path = os.path.join(self.path, fn)
                im = PostScriptImage(complete_path, self.out_dir, self.path)
                self.images.append(im)
                self.images_by_name[base_fn] = im
    def get_sorted_images(self):
        keys = self.images_by_name.keys()
        keys.sort()
        return [self.images_by_name[key] for key in keys]
    def get_html(self, out_dir):
        tmpl = open("template/results-dir.html", "r").read()
        images = "<ul class=\"results_charts\">\n"
        for image in self.get_sorted_images():
            images += "<li>" + image.get_html() + "</li>\n"
        images += "</ul>\n"
        files_html = ""
        for txt_file in self.txt_tables:
            t = TextTable()
            t.read_file(os.path.join(self.path, txt_file))
            files_html += """
            <h5>%s</h5>
            %s
            """ % (txt_file, t.get_html())
        files_html += "<ul>\n"
        for t in self.txt_files:
            files_html += "<li>%s</li>\n" % t.get_html()
        files_html += "</ul>\n"
        return tmpl % {'title': self.get_pretty_name(), 
            'images': images, 
            'files': files_html,
            'anchor': self.get_slug(),
            'date': self.get_execution_time(), }
    def get_execution_time(self):
        t = self.log.get_execution_time()
        if not t:
            t = "Completion time missing"
        return t
    def get_slug(self):
        return slugify(self.path)
    def copy_files(self):
        for f in self.txt_files:
            shutil.copy(f.get_src_fn(), f.get_dst_fn())
            os.chmod(f.get_dst_fn(), stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH | stat.S_IWUSR)
    def get_order_key(self):
        r = re.search(r'chr(?P<foo>[0-9]+)', self.path)
        if not r:
            return self.path
        else:
            return int(r.groupdict()['foo'])
    def get_pretty_name(self):
        return "Chromosome %s" % (self.get_order_key(), )


def main():
    parser = OptionParser()
    parser.add_option("-o", "--output", dest="dest_dir",
            help="write report to DIR", metavar="DIR")
    parser.add_option("-d", "--results-parent-dir", dest="results_parent_dir",
            help="parent directory for results", metavar="DIR")
    (options, args) = parser.parse_args()

    assert options.dest_dir, "Please specify destination directory"
    assert options.results_parent_dir, "Please results parent directory"

    rds = ResultsDirContainer(options.dest_dir)

    for results_dir in os.listdir(options.results_parent_dir):
        results_path = os.path.join(options.results_parent_dir, results_dir)
        if os.path.isdir(results_path):
            rds.add_dir(results_path)

    rds.generate_html()

if __name__ == '__main__':
    main()


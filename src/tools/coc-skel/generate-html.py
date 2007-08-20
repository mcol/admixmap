#!/usr/bin/python
# -*- coding: UTF-8 -*-
# 
# HTML generator for Hapmixmap results
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
import os, os.path
from optparse import OptionParser

from hapmixmap import *

# Elements to display:
#
# - Coefficient of constraint
# - Posterior Quantiles
# - Energy plots (all of them?)
# - Arguments file

def main():
    parser = OptionParser()
    parser.add_option("-o", "--output", dest="dest_dir", help="write report to DIR", metavar="DIR")
    parser.add_option("-c", "--config", dest="config", help="configuration file", metavar="FILE")
    # parser.add_option("-d", "--results-parent-dir", dest="results_parent_dir", help="parent directory for results", metavar="DIR")
    parser.add_option("-m", "--main-dir", dest="main_dir", help="Main directory", metavar="DIR")
    (options, args) = parser.parse_args()

    assert options.main_dir, "Please specify main directory"
    assert options.dest_dir, "Please specify destination directory"
    assert options.config, "Please specify configuration file"


    main_dir = os.path.join("..", options.main_dir)
    results_group = AnalysisResultsGroup(main_dir, options.dest_dir)
    for config_line in open(options.config, "r").readlines():
        if config_line[0] == '#':
            continue
        config_tokens = config_line.split(':')
        src_dir = config_tokens[0]
        results_group.add_src_directory(src_dir)

    results_group.write_html()
    results_group.write_csv()

if __name__ == '__main__':
    main()


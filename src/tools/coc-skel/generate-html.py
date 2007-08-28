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


#!/usr/bin/python
# -*- coding: utf-8 -*-
# 
# Submit a taskfarm job to Walton
#
# Maciej Bliziński

import sys
import os
import os.path
import re
from optparse import OptionParser

usage = False
parser = OptionParser("usage: %prog options")

parser.add_option("-f", "--file", dest = "file_name",
        help = "Script to submit", metavar = "FILE")
parser.add_option("-t", "--time", dest = "hours",
        help = "Time to run in hours")
parser.add_option("-n", "--nodes", dest = "nodes",
        help = "Number of nodes (2 procs per node)")
parser.add_option("-e", "--template", dest = "template",
        help = "PBS script template")
parser.add_option("-s", "--submit", dest = "submit",
        action = "store_true",
        help = "Submit the job to PBS")

options, args = parser.parse_args()

if not options.template:
    print "Please specify template file name."
    usage = True
if not options.file_name:
    print "Please specify script file name."
    usage = True
if not options.hours:
    print "Please specify number of hours."
    usage = True
if not options.nodes:
    print "Please specify number of nodes."
    usage = True

if usage:
    print parser.format_help()
    sys.exit(1)

# template = open("taskfarm-template.pbs").read()
template = open(options.template).read()

work_dir = os.path.abspath(".")

if re.search(r':', options.hours):
    hours, minutes = options.hours.split(':')
else:
    hours = options.hours
    minutes = '00'

pbs_script = template % {
    'name': options.file_name,
    'nodes': options.nodes,
    'hours': hours,
    'minutes': minutes,
    'script': options.file_name,
    'workdir': work_dir, }

open("tmp.pbs", "w").write(pbs_script)

if options.submit:
    os.system("qsub tmp.pbs")
else:
    print "PBS Script would be:"
    print
    print pbs_script

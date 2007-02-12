#!/usr/bin/python
# 

import re
import os.path
from urllib import FancyURLopener

from hapmap import HapmapChromosomeLegend

# 1. Download all the legend files
# 2. Build a dictionary of rs ids
# 3. Save the dictionary on the disk.
# 4. Given a list of loci, return a chromosome URL or return an error.

BASE_URL = 'http://www.hapmap.org/downloads/phasing/2006-07_phaseII/phased/'
WORKING_DIR = 'legend'

def download_legends(base_url, working_dir):
    print "Retrieving legend files list."
    print BASE_URL
    uo = FancyURLopener()
    file_list = uo.open(base_url)
    lines = file_list.readlines()

    legend_lines = []
    for line in lines:
        tokens = re.split(r'"', line)
        legend_lines.extend(filter(lambda(x): re.search(r'legend\.txt', x), tokens))

    legend_lines = filter(lambda(x): not(re.search(r'>', x)), legend_lines)
    legend_lines = list(set(legend_lines))
    legend_lines.sort()
    # Download all the legend files
    counter = 0L
    for legend_file_name in legend_lines:
        local_file_name = "%s/%s" % (working_dir, legend_file_name)
        if not os.path.exists(local_file_name):
            print "Downloading %s (%s of %s)" % (legend_file_name, counter, len(legend_lines))
            lf = uo.open("%s%s" % (base_url, legend_file_name))
            open(local_file_name, "w").write(lf.read())
        counter += 1
    print "All files downloaded."

download_legends(BASE_URL, WORKING_DIR)

# List the files in the directory.
for root, dirs, files in os.walk(WORKING_DIR):
    # print root, dirs, files
    files.sort()
    for file_name in files:
        pass
        # lf = HapmapChromosomeLegend(os.path.join(WORKING_DIR, file_name))
        # lf.test()




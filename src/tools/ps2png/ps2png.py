#!/usr/bin/python
# -*- coding: UTF-8 -*-
# 
# ps2png.py, Postscript to PNG converter.
#
# Requirements:
# imagemagick
# ghostscript
#
# This script is a wrapper around imagemagick, which uses ghostscript
# to read postscript.

import os
import os.path
import sys

def convert_ps_to_png(file_name):
    new_file_name = file_name[:-3] + ".png"
    print "convert_ps_to_png(%s) -> '%s'" % (file_name, new_file_name)
    os.system("convert \"%s\" -rotate 90 \"%s\"" % (file_name, new_file_name))

def main():
    for root, dirs, files in os.walk('.'):
        for file in files:
            if file[-3:] == '.ps':
                convert_ps_to_png(file)

if __name__ == '__main__':
    main()

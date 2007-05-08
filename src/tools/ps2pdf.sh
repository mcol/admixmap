#!/bin/bash

gs -dNOPAUSE -sDEVICE=pdfwrite -dTextAlphaBits=4 -dGraphicsAlphaBits=4 -r110x110 -sOutputFile=$1.pdf -q $1.ps -c quit 
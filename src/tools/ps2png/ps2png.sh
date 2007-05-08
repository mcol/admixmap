#!/bin/bash

gs -dNOPAUSE -sDEVICE=png256 -dTextAlphaBits=4 -dGraphicsAlphaBits=4 -r110x110 -sOutputFile=$1.png -q $1.ps -c quit 
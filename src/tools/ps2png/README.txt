ps2png.py

A Python script which converts all the Postscript files from the current
directory into PNG format.

Requirements:
 * Python
 * Imagemagick
 * Ghostscript

Direct links to download Windows binaries:
 * http://www.python.org/ftp/python/2.5.1/python-2.5.1.msi
 * http://www.imagemagick.org/download/binaries/ImageMagick-6.3.3-10-Q16-windows-static.exe
 * http://mirror.switch.ch/ftp/mirror/ghost/GPL/gs856/gs856w32.exe


ps2png.sh
Uses ghostscript to convert a Postscript file to PNG format.
Usage: ./ps2png filename
NB: omit the .ps extension

ps2png.bat
DOS equivalent of the above. Assumes Ghostscript is installed in C:\gs\gs8.11
NB: omit .ps extensions
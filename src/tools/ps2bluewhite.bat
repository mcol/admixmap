:: DOS batch script to change black-on-white
:: postscript plots to white-on blue

:: change background colour from unspecified to blue
perl -p -i.bak -e "s/bp\n/bp\n\/bg { 0 0 1 } def\n0.00 0.00 841.89 595.28 r p2\n/g;" %1

:: change colours of plotting elements (points, axes etc) from black to white
perl -p -i.bak -e "s/0 0 0 rgb/1 1 1 rgb/g;" %1
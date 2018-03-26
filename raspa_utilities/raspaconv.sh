#!/bin/bash

egrep "absolute adsorption. ......... .avg" Output/*/* > aa

cat >> aa.gnuplot << EOF
set xlabel "printed step"
set ylabel "mol/uc"
  p "aa" u 3 w l  lc rgb 'red'   ti "current step value"
rep "aa" u 5 w l  lc rgb 'blue'  ti "average value"
pause -1
EOF

gnuplot aa.gnuplot 

rm -rf aa aa.gnuplot

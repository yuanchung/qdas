#!/bin/sh

jname=$1

gnuplot << EOF

set style line 1 lw 3 lt 1 pt 6    # solid line, circle point
set style line 2 lw 3 lt 11 pt 4   # long dash line, square point
set style line 3 lw 3 lt 5 pt 12   # dash-dot line, diamond point
set style line 4 lw 3 lt 0 pt 3    # dot line, cross point
set style line 5 lw 3 lt 2 pt 8    # thin-long dash line, triangular point
set style line 6 lw 3 lt 3 pt 14   # short dash line,  pentagon point
set style line 7 lw 3 lt 23 pt 2    # long dash-dot line, cross point

set terminal postscript portrait enhanced color

set encoding iso_8859_1
set output "${jname}.ps"

set key off
set multiplot
set size 1,0.5

set origin 0,0
set logscale x
set title "Model RCOX 3PPE peak shift"
set ylabel "peak shift (fs)"
set xlabel "T (fs)"
plot [1:1000] [] 'rcox_peps.dat' u 1:2 t ""  w l

set origin 0,0.5
unset logscale x
set logscale y
set title "Model RCOX 3PPE signal, (Set ${jname}, 2C downhill)"
set xlabel "{/Symbol t} (fs)"
set ylabel "T (fs)"
set pm3d
set isosample 40
set contour base
set cntrparam levels 10
set view 0,0
splot [] [10:1000] 'rcox.ndat' u 1:2:3 t ""  w l

set nomultiplot

EOF


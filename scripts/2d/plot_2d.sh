#!/bin/sh

# plot real part 2D spectrum

# Usage:
# ./plot_2d.sh three_0325-2d.dat

for i in $*
do

echo "Processing $i ..."

# assume the filename is in the prefix_T-2d.dat form, we get the pop. time T
T=`echo $i | cut -d'-' -f1 | tr -d 'a-zA-Z_'`
T=`expr ${T} + 1 - 1` # this get rid of leading zeros

dname=`basename $PWD`
bname=`basename ${i} .dat`
bname="${dname}_${bname}"

# call gnuplot to plot the 2D map
gnuplot << EOF

set terminal postscript eps enhanced color solid 16
set encoding iso_8859_1
set output "${bname}_map.eps"

set label 1 "{/Symbol w}_{/Symbol t} (cm^{-1})" at screen 0.47,0.1
set label 2 "{/Symbol w}_t (cm^{-1})" at screen 0.19,0.45 rotate by 90
set label 3 "T=${T} fs" at screen 0.6,0.78 font "Helvetica,20" front

set pm3d map
set size square
set xrange [50:700]
set yrange [50:700]
#set zrange [-6:6]

set palette rgbformulae 22,13,-31
set cbrange [-4:4]

set contour surface
set style data lines
set cntrparam levels discrete 7,6,5,4,3,2,1,0.5,-0.5,-1,-2,-3,-4,-5,-6,-7

unset key
unset colorbox
unset xlabel
unset ylabel
unset clabel

set multiplot


set pm3d map; set surface
splot '${i}' u 1:2:(\$3/1e6)
unset surface
splot '${i}' u 1:2:(\$3/1e6) w l lw 0.3 lc rgb "#666666"

unset multiplot

EOF

eps2eps ${bname}_map.eps aaa.eps
mv aaa.eps ${bname}_map.eps

done


#!/bin/bash

# plot slices of 2D spectrum

PROG_2DANALYSIS=~/x-prog/qdas/tools/2danalysis

# standard file; the 2d plots will be scaled to the maximum of
# the real part of this spec.
STANDARD=dimer_0000-2d.dat

# define how to plot diagonal and vertical slices in an array; format is
# filename_prefix:xlabel:ylabel:datalabel ...
VSLICE_ARRAY[1]="2d-d:T (fs):{/Symbol w}_t/{/Symbol w}_{/Symbol t} (cm^{-1}):diag. (total)"
VSLICE_ARRAY[2]="2d-ptau-d:T (fs):{/Symbol w}_t/{/Symbol w}_{/Symbol t} (cm^{-1}):diag. (rephasing)"
VSLICE_ARRAY[3]="2d-ntau-d:T (fs):{/Symbol w}_t/{/Symbol w}_{/Symbol t} (cm^{-1}):diag. (non-rephasing)"
VSLICE_ARRAY[4]="2d-v250:T (fs):{/Symbol w}_t (cm^{-1}) (cm^{-1}):{/Symbol w}_{/Symbol t}=250 (total)"
VSLICE_ARRAY[5]="2d-ptau-v250:T (fs):{/Symbol w}_t (cm^{-1}):{/Symbol w}_{/Symbol t}=250 (rephasing)"
VSLICE_ARRAY[6]="2d-ntau-v250:T (fs):{/Symbol w}_t (cm^{-1}):{/Symbol w}_{/Symbol t}=250 (non-rephasing)"

# define how to plot horizontal slices in an array; format is
# filename_prefix:xlabel:ylabel:datalabel ...
HSLICE_ARRAY[1]="2d-h250:{/Symbol w}_{/Symbol t} (cm^{-1}):T (fs):{/Symbol w}_t=250 (total)"
HSLICE_ARRAY[2]="2d-ptau-h250:{/Symbol w}_{/Symbol t} (cm^{-1}):T (fs):{/Symbol w}_t=250 (rephasing)"
HSLICE_ARRAY[3]="2d-ntau-h250:{/Symbol w}_{/Symbol t} (cm^{-1}):T (fs):{/Symbol w}_t=250 (non-rephasing)"

# nothing to be changed after this

### no need to change stuff after this line

# scaling factor for graphs; data points will be divided by this factor
# before the plot
SCALE=`$PROG_2DANALYSIS rmax ${STANDARD} | grep 'Maximum real part' | cut -d'=' -f2 | tr -d ' '`

echo "Will use ${STANDARD} as the standard and scale ${SCALE} to 1.0"
echo " "

# a function that plots a 2d map for diag. and vertical slices
# usage plot_vslice_data data_file output_file xlabel ylabel data_label
function plot_vslice_data {

# get arguments
    datafile=$1
    outfile=$2
    xlabel=$3
    ylabel=$4
    dlabel=$5

# call gnuplot to plot the 2D map
    gnuplot << EOF
set terminal postscript eps enhanced color solid 16
set encoding iso_8859_1
set output "$outfile"

set label 1 "$xlabel" at screen 0.47,0.1
set label 2 "$ylabel" at screen 0.19,0.45 rotate by 90
set label 3 "$dlabel" at screen 0.72,0.78 right font "Helvetica,20" front

set pm3d map
set size square
set xrange [0:600]
set yrange [0:1000]
#set zrange [-6:6]

set palette rgbformulae 22,13,-31
set cbrange [-0.9:0.9]

set contour surface
set style data lines
set cntrparam levels discrete 1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8,-0.9,-1

unset key
unset colorbox
unset xlabel
unset ylabel
unset clabel

set multiplot


set pm3d map; set surface
splot '$datafile' u 1:2:(\$3/$SCALE)
unset surface
splot '$datafile' u 1:2:(\$3/$SCALE) w l lw 0.1 lc rgb "#666666"

unset multiplot

EOF

    eps2eps $outfile $$.aaa.eps
    mv $$.aaa.eps $outfile
} # end of function plot_vslice_data

# a function that plots a 2d map for horizontal slices
# usage plot_hslice_data data_file output_file xlabel ylabel data_label
function plot_hslice_data {

# get arguments
    datafile=$1
    outfile=$2
    xlabel=$3
    ylabel=$4
    dlabel=$5

# call gnuplot to plot the 2D map
    gnuplot << EOF
set terminal postscript eps enhanced color solid 16
set encoding iso_8859_1
set output "$outfile"

set label 1 "$ylabel" at screen 0.19,0.45 rotate by 90
set label 2 "$xlabel" at screen 0.47,0.1
set label 3 "$dlabel" at screen 0.72,0.79 right font "Helvetica,16"

set pm3d map
set size square
set yrange [600:0]
set xrange [0:1000]
#set zrange [-6:6]

set palette rgbformulae 22,13,-31
set cbrange [-0.9:0.9]

set contour surface
set style data lines
set cntrparam levels discrete 1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8,-0.9,-1

unset key
unset colorbox
unset xlabel
unset ylabel
unset clabel

set multiplot


set pm3d map; set surface
splot '$datafile' u 2:1:(\$3/$SCALE)
unset surface
splot '$datafile' u 2:1:(\$3/$SCALE) w l lw 0.1 lc rgb "#666666"

unset multiplot

EOF

    eps2eps $outfile $$.aaa.eps
    mv $$.aaa.eps $outfile
} # end of function plot_hslice_data


## main program

# plot vertical and diagonal slices
for i in `seq 1 ${#VSLICE_ARRAY[*]}`
do
    pstr=${VSLICE_ARRAY[$i]}
    prefix=`echo $pstr | cut -d':' -f1`
    label1=`echo $pstr | cut -d':' -f2`
    label2=`echo $pstr | cut -d':' -f3`
    label3=`echo $pstr | cut -d':' -f4`
    echo "Plotting ${prefix}.dat via \"$pstr\" ..."
    plot_vslice_data ${prefix}.dat ${prefix}.eps "${label1}" "${label2}" "${label3}"
done

# plot horizontal slices
for i in `seq 1 ${#HSLICE_ARRAY[*]}`
do
    pstr=${HSLICE_ARRAY[$i]}
    prefix=`echo $pstr | cut -d':' -f1`
    label1=`echo $pstr | cut -d':' -f2`
    label2=`echo $pstr | cut -d':' -f3`
    label3=`echo $pstr | cut -d':' -f4`
    echo "Plotting ${prefix}.dat via \"$pstr\" ..."
    plot_hslice_data ${prefix}.dat ${prefix}.eps "${label1}" "${label2}" "${label3}"
done

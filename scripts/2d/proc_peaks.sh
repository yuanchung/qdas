#!/bin/sh

# process all *-2d.dat files and use 2danalysis to find integrated volume 
# under cross peaks

PROG_2DANALYSIS=~/x-prog/qdas/tools/2danalysis

# define regions of peaks, the format is
# filename_prefix,w1lower,w1upper,w2lower,w2upper
PEAKLIST="peak21,290,420,150,250 peak31,450,650,150,250 peak32,450,650,290,420 peak12,150,250,290,420 peak13,150,250,450,650 peak23,290,420,450,650 peak11,150,250,150,250 peak22,290,420,290,420 peak33,450,650,450,650"

# nothing to be changed after this

for i in $PEAKLIST
do
  prefix=`echo $i | cut -d',' -f1`
  echo "Processing $prefix ..."
  area=`echo $i | cut -d',' -f2- | tr ',' ' '`
  echo "  peak area defined by ($area)."
  rm -f  ${prefix}.out
  touch ${prefix}.out
  echo "  processing ans save 2danalysis output in ${prefix}.out."
  for j in *-2d.dat
  do
    $PROG_2DANALYSIS $j integrate $area >> ${prefix}.out
  done
  echo "  save data points in ${prefix}.dat."
  grep 'Integrated volume' ${prefix}.out | \
       sed 's/-2d.dat//g' | tr -s 'a-zA-Z_():=,' ' ' > ${prefix}.dat
done



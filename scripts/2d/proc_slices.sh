#!/bin/sh
#
# use 2danalysis to find slices of 2D spectra (for all T)
#

PROG_2DANALYSIS=~/x-prog/qdas/tools/2danalysis

# define w1 for vertical slices, the format is
# filename_prefix,w1
VSLICE_LIST="v250,250"

# define w2 for horizontal slices, the format is
# filename_prefix,w2
HSLICE_LIST="h250,250"

# nothing to be changed after this

# we loop over *-ntau-2d.dat, *-ptau-2d.dat, ...
SLIST="2d.dat ntau-2d.dat ptau-2d.dat"

for postfix in $SLIST
do

FLIST=`ls *[0-9]-${postfix}`
FPREFIX=`basename 2d-${postfix} -2d.dat`

## get the diagonal slices
echo "Processing and saving diagonal slice output in ${FPREFIX}-d.out."
rm -f ${FPREFIX}-d.out
touch ${FPREFIX}-d.out
for j in $FLIST
do
  $PROG_2DANALYSIS dslice $j >> ${FPREFIX}-d.out
done
echo "  saving diagonal slice into ${FPREFIX}-d.dat ..."
grep "Diagonal slice" ${FPREFIX}-d.out | \
       sed 's/-2d.dat//g' | tr -s 'a-zA-Z_():=,' ' ' | cut -d' ' -f2,3,5,6 \
       > ${FPREFIX}-d.dat

# vertical slices
for i in $VSLICE_LIST
do
  prefix=`echo $i | cut -d',' -f1`
  echo "Processing $prefix and saving 2danalysis output in ${FPREFIX}-${prefix}.out ..."
  omega=`echo $i | cut -d',' -f2 | tr ',' ' '`
  echo "  slice defined by w1=$omega."
  rm -f  ${FPREFIX}-${prefix}.out
  touch ${FPREFIX}-${prefix}.out
  for j in $FLIST
  do
    $PROG_2DANALYSIS vslice $omega $j >> ${FPREFIX}-${prefix}.out
  done
  echo "  save data points in ${FPREFIX}-${prefix}.dat"
  grep 'Vertical slice' ${FPREFIX}-${prefix}.out | \
       sed 's/-2d.dat//g' | tr -s 'a-zA-Z_():=,' ' ' | cut -d' ' -f2,4,5,6 \
       > ${FPREFIX}-${prefix}.dat
done

# horizontal slices
for i in $HSLICE_LIST
do
  prefix=`echo $i | cut -d',' -f1`
  echo "Processing $prefix and saving 2danalysis output in ${FPREFIX}-${prefix}.out ..."
  omega=`echo $i | cut -d',' -f2 | tr ',' ' '`
  echo "  slice defined by w2=$omega."
  rm -f  ${FPREFIX}-${prefix}.out
  touch ${FPREFIX}-${prefix}.out
  for j in $FLIST
  do
    $PROG_2DANALYSIS hslice $omega $j >> ${FPREFIX}-${prefix}.out
  done
  echo "  save data points in ${FPREFIX}-${prefix}.dat"
  grep 'Horizontal slice' ${FPREFIX}-${prefix}.out | \
       sed 's/-2d.dat//g' | tr -s 'a-zA-Z_():=,' ' '| cut -d' ' -f2,3,5,6 \
       > ${FPREFIX}-${prefix}.dat
done

done # done the postfix loop

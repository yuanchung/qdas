#!/bin/sh
#
# Script to compute peakshift results from data produced by
# jobs created from gen_2d.sh. Note that for simlicity we 
# only process positive tau data for peak-shifts.
#

# we need cspline_max and normalize; those are in the ../tools dir
CSPLINE_MAX=~/tools/cspline_max
GAU_MAX=~/tools/gau_max
NORMALIZE=~/tools/normalize

# A prefix is required as the argument,
# output files is named as prefix_tau_T.out

if [ "$#" != "1" ] ; then
  echo "Filename prefix of .out files is required!!"
  exit
else
  prefix=$1
  echo "Processing .out files named as ${prefix}_tau_T.out..."
fi

str=`find . -name "${prefix}_[0-9]*_*.out" 2>/dev/null | head`
if [ -z "$str" ] ; then
  echo "Output files with prefix \"${prefix}\" not found!!"
  exit
fi

# backup old data file to make space for new files
echo "Will save integrated signals in ${prefix}.dat"
if [ -e ${prefix}.dat ] ; then
  echo "  Found old .dat file, backup to ${prefix}.dat-bak"
  mv ${prefix}.dat ${prefix}.dat-bak
fi
echo "Will save normalized integrated signals in ${prefix}.ndat"
if [ -e ${prefix}.ndat ] ; then
  echo "  Found old .ndat file, backup to ${prefix}.ndat-bak"
  mv ${prefix}.ndat ${prefix}.ndat-bak
fi
echo "Will save peakshift data in ${prefix}_peps.dat"
if [ -e ${prefix}_peps.dat ] ; then
  echo "  Found old _peps.dat file, backup to ${prefix}_peps.dat-bak"
  mv ${prefix}_peps.dat ${prefix}_peps.dat-bak
fi

# gather list of population times
TLIST=`find . -name "${prefix}_[0-9]*_*.out" | cut -d'_' -f3 | sed 's/.out//g' | sort -n | uniq`

# note that the definition of population time in negative tau 2d jobs 
# are different from that in the PEPS jobs. So we first make links and 
# adjust population time on the file names.
echo "Transforming population time for negative coherence time jobs..."
for i in $TLIST
do
  # gather a list of negative 2D tau for this T set
  NTAULIST=`ls ${prefix}_-[0-9]*_${i}.out | cut -d'_' -f2 | sort -n | uniq`
  for j in $NTAULIST
  do
    # PEPS T' = T-tau
    new_T=`expr ${i} - ${j}`
    new_T_str=`printf "%04d" ${new_T}`
    ln -sf ${prefix}_${j}_${i}.out ${prefix}_${j}_${new_T_str}.out-tmp$$
  done
done

# now we are ready to process each T job set
for i in $TLIST
do
  echo "Processing T=${i} data..."
  grep 'Int. Signal' ${prefix}_[0-9]*_${i}.out \
                           ${prefix}_-[0-9]*_${i}.out-tmp$$ 2>/dev/null | \
       tr '=' ',' | cut -d',' -f2,4,6 | tr -d ',' | sort -n > $$.tmp
  cat $$.tmp >> ${prefix}.dat # pile up data in ${prefix}.dat
  echo >> ${prefix}.dat
  # split the data in two column, normalize it, then save
  cat $$.tmp | tr -s ' ' ' ' | cut -d' ' -f2,4 | $NORMALIZE > $$-norm.tmp
  cat $$.tmp | sed 's/[-0-9.]*$//g' > $$-1.tmp # tau and T pairs
  cat $$-norm.tmp | tr -s ' ' ' ' | cut -d' ' -f3 > $$-2.tmp
  paste $$-1.tmp $$-2.tmp >> ${prefix}.ndat
  echo >> ${prefix}.ndat
  # now the peakshift data
  peps1=`cat $$-norm.tmp | $CSPLINE_MAX | cut -d' ' -f1`
  peps2=`cat $$-norm.tmp | $GAU_MAX | cut -d' ' -f1`
  echo "${i}  ${peps1} ${peps2}" >> ${prefix}_peps.dat
  rm -f $$.tmp $$-norm.tmp $$-1.tmp $$-2.tmp
done

# remove tmp links
find . -maxdepth 1 -name "*.out-tmp$$" -exec rm {} \;


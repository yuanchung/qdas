#!/bin/sh

# process tnl-dm3pes output and generate 2D spectra.

# A prefix is required as the argument,
# and tnl-dm3pes output files is named as prefix_tau_T.out

# we need 2dspec; those are in the ../tools dir
PROG_2DSPEC=~/tools/2dspec

# A prefix is required as the argument,
# output files is named as prefix_tau_T.out

if [ "$#" != "1" ] ; then
  echo "Filename prefix of .out files is required!!"
  exit
else
  prefix=$1
  echo "Processing .out files named as ${prefix}_tau_T.out..."
fi

str=`find . -name "${prefix}_*_*.out" 2>/dev/null | head`
if [ -z "$str" ] ; then
  echo "Output files with prefix \"${prefix}\" not found!!"
  exit
fi

# backup old data file to make space for new files
echo "Will save P(tau,t) in ${prefix}_T.tt ..."
echo "Will save 2d spectra in ${prefix}_T.2d ..."

# gather list of population times
TLIST=`find . -name "${prefix}_*_*.out" | cut -d'_' -f3 | sed 's/.out//g' | sort -n | uniq`

for T in $TLIST
do
  echo "Processing T=${T} data..."

  # gather list of coherence times
  TAULIST=`ls ${prefix}_*_${T}.out | cut -d'_' -f2 | sort -n | uniq`
  NTAU=`ls ${prefix}_*_${T}.out | cut -d'_' -f2 | sort -n | uniq | wc -l`
  # collect P(tau,t)
  rm -f ${prefix}_${T}.tt
  for tau in $TAULIST
  do
    Nt=`grep -A 4000 -E "^t=[- ]+0\.0000," ${prefix}_${tau}_${T}.out | \
        grep P= | wc -l`
    grep -A 4000 -E "^t=[- ]+0\.0000," ${prefix}_${tau}_${T}.out | \
        grep P= | sed "s/^t=/${tau}/g" | tr -d 'a-zA-Z=+|,' | tr -s ' ' ' ' | \
	cut -d' ' -f1,2,3,4 >> ${prefix}_${T}.tt
    echo " " >> ${prefix}_${T}.tt # this is for gnuplot's splot
  done
  # now ready to do FFT
  #echo "  2D FFT on ${prefix}_${T}.tt ..."
  #$PROG_2DSPEC ${NTAU} ${Nt} ${prefix}_${T}.tt > ${prefix}_${T}.2d
done


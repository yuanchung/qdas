#!/bin/sh

# use 2dspec to convert .tt file to 2d spectra
PROG_2DSPEC="2dspec"

# number of tau and t points in the input file
if [ $# -ge 3 ] ; then
  N_TAU=$1
  shift
  N_t=$1
  shift
else
  N_TAU=161
  N_t=81
fi

# size of fft
FFT_SIZE=512

# Window for output
W1_RANGE="0 1200"
W2_RANGE="0 1200"

# number of positive/negative tau points; assume symmetry
HALF_N_TAU=`expr ${N_TAU} / 2 + 1`

for i in $*
do 
  echo "Processing $i ..."
  bname=`basename $i .tt`
  $PROG_2DSPEC -w1range $W1_RANGE -w2range $W2_RANGE \
               -fftsize ${FFT_SIZE} ${FFT_SIZE} \
               ${N_TAU} ${N_t} $i > ${bname}-2d.out
  grep Sr= ${bname}-2d.out | tr -d "a-zA-Z=+:;,][|\\/_" > ${bname}-2d.dat
  #grep Sr_norm= ${bname}-2d.out | tr -d "a-zA-Z=+:;,][|\\/_" > ${bname}-2d.ndat

  # uncomment the following part if separately handeling the 
  # positive tau and negative tau parts are required.
  
  echo "  rephasing part in ${bname}-ptau.tt ..."
  grep -A 1000000 '^000 -*0.0000' $i > ${bname}-ptau.tt
  $PROG_2DSPEC -w1range $W1_RANGE -w2range $W2_RANGE \
               -fftsize ${FFT_SIZE} ${FFT_SIZE} \
               ${HALF_N_TAU} ${N_t} ${bname}-ptau.tt > ${bname}-ptau-2d.out
  grep Sr= ${bname}-ptau-2d.out | tr -d "a-zA-Z=+:;,][|\\/_" \
		> ${bname}-ptau-2d.dat
  #grep Sr_norm= ${bname}-ptau-2d.out | tr -d "a-zA-Z=+:;,][|\\/_" \
  #		> ${bname}-ptau-2d.ndat

  echo "  non-rephasing part in ${bname}-ntau.tt ..."
  last_t=`tail $i | cut -d' ' -f2 | sort -g | tail -n 1`
  grep -B 1000000 "^000 ${last_t}" $i > ${bname}-ntau.tt
  echo "" >> ${bname}-ntau.tt
  $PROG_2DSPEC -w1range $W1_RANGE -w2range $W2_RANGE \
               -fftsize ${FFT_SIZE} ${FFT_SIZE} \
               ${HALF_N_TAU} ${N_t} ${bname}-ntau.tt > ${bname}-ntau-2d.out
  grep Sr= ${bname}-ntau-2d.out | tr -d "a-zA-Z=+:;,][|\\/_" \
		> ${bname}-ntau-2d.dat
  #grep Sr_norm= ${bname}-ntau-2d.out | tr -d "a-zA-Z=+:;,][|\\/_" \
  #		> ${bname}-ntau-2d.ndat

  # remove .out files...
  rm -f ${bname}-2d.out ${bname}-ptau-2d.out ${bname}-ntau-2d.out

done


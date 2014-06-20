#!/bin/sh

# extract the negative tau part (non-rephasing part) of a tau-t map
# from a full .tt file; used to generate the non-rephasing-only
# 2d spectrum

last_t=`tail $1 | cut -d' ' -f2 | sort -n | tail -n 1`

grep -B 1000000 "^000 ${last_t}" $1 > $$.tmp
echo "" >> $$.tmp

cat $$.tmp

rm $$.tmp



#!/bin/sh

# extract the positive tau part (rephasing part) of a tau-t map
# from a full .tt file; used to generate the rephasing-only
# 2d spectrum

grep -A 1000000 '^000 -*0.0000' $1



#!/bin/bash

# qs_archive.sh PREFIX

# make an archive of all .done .out files

if [ "z$1" = "z" ] ; then
  echo "Usage: qs_archive.sh PREFIX"
  exit
fi

echo "Processing $1_archive.tar.gz ..."
tar cfz $1_archive.tar.gz $1_*.key.done $1_*.out

echo "Removing $1_*.key.done and $1_*.out files..."
rm -f $1_*.key.done $1_*.out


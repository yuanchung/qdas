#!/bin/sh
#
# Update all 2cdh 3PEPS plots
#

PLOT_RC="${PWD}/plot_rcox.sh"
FLIST=`find . -name rcox_peps.dat`
basedir=`pwd`

for i in $FLIST
do
  dirname=`dirname ${i}`
  jname=`dirname ${i} | sed 's/\.\///' | tr '/_' '-'`
  echo "Plotting $jname ..."
  cd ${dirname}
  $PLOT_RC ${jname}
  cd $basedir
done



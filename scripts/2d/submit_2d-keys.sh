#!/bin/sh

# Script that submits LSF jobs for a 2D spectra run

# resource and queue setup
queue=short

jname=`basename $PWD`
errfile=${jname}.err
logfile=${jname}.log

email=yccheng@berkeley.edu

# submit a job for each keyword files

for i in `ls -t *.key`
do
  bname=`basename ${i} .key`
  if [ ! -f ${bname}.out ] ; then
  bsub -B -N -u ${email} -q ${queue} -J ${jname}-${i} -o ${logfile} -e ${errfile} $HOME/tools/tnl-dm3pes_single.sh $i
  else
  echo "Skip $i ..."
  fi
  sleep 3
done


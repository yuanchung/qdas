#!/bin/sh
#
# Generate RC 3PPEPS jobs from a tmpl file, this script with correct
# supporting files can be used to setup a series of 3PPES tnl-dm3pes
# Jobs.

# required files:
# ${JOBNAME}.key-tmpl-tmpl: a template of template tnl-dm3pes keyword file
# ${JOBNAME}.pbs-tmpl: template PBS job file.
# gen_3ppes.sh: script use to generate .key and .pbs files.

# usage: to use this script, modify JOBNAME and gamma lists below.
#        also modify JOBNAME and 3ppes scan ranges in the gen_3ppes.sh script.
# ./run_3ppes_jobs.sh

JOBNAME="rcox"
JOBDIRNAME="rcox_peps"

# which PBS queue to submit to?
PBS_QUEUE="short"

# gamma lists, these number will be sub. into ${JOBNAME}.key-tmpl-tmpl

GBL_LIST="0.5 1.0 1.5 2.0"
GHL_LIST="0.5 1.0 1.5 2.0"
GHB_LIST="0.5 1.0 2.0"

# no parameter to be changed below

WORK_DIR=`pwd`

for gbl in $GBL_LIST
do
  for ghl in $GHL_LIST
  do
    for ghb in $GHB_LIST
    do
	  dname="${JOBDIRNAME}-${gbl}-${ghl}-${ghb}"
	  if [ -d "$dname" ] ; then
	    echo "Error!! $dname already exists!!"
	    echo "  Skip $dname!!"
	    continue
	  fi
	  echo "Preparing $dname ..."
	  mkdir $dname
          cd $dname
	  cp ../gen_3ppes.sh ../${JOBNAME}.pbs-tmpl .
          cat ../${JOBNAME}.key-tmpl-tmpl | \
	      sed "s/VAR_GBL/ $gbl /g" | \
	      sed "s/VAR_GHL/ $ghl /g" | \
	      sed "s/VAR_GHB/ $ghb /g" > ${JOBNAME}.key-tmpl
	  # generate .key files and then submit the job
	  echo "  Generating keyword files..."
	  ./gen_3ppes.sh > /dev/null
	  echo "  Submitting the job to PBS..."
          qsub -q $PBS_QUEUE *.pbs
	  
          # This job down, do the next one
	  cd ..
    done
  done
done

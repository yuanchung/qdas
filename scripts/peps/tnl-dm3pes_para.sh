#!/bin/sh
#
# A wrapper script used to run parallel tnl-dm3pes jobs
#
# Usage: ssh nodename tnl-dm3pes_para.sh nodefile_string work_dir keyfile
#

# Global setup
TNL_DM3PES=/home/platin/tools/tnl-dm3pes

node=$1
workdir=$2
keyfile=$3

jname=`basename ${keyfile} .key`

# run job
cd $workdir
${TNL_DM3PES} ${keyfile} > ${jname}.out
# mark the job done, this makes the PBS job re-runable
mv ${jname}.key ${jname}.key.done

# Job done, set the node as a free node
hostname > ${node}.$$.freenode


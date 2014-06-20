#!/bin/sh
#
# A wrapper script used to run a single tnl-dm3pes job
#
# Usage: tnl-dm3pes_single.sh keyfile
#

# Global setup
TNL_DM3PES=$HOME/tools/tnl-dm3pes

keyfile=$1

jname=`basename ${keyfile} .key`

# run job
${TNL_DM3PES} ${keyfile} > ${jname}.out
# mark the job done, this makes the PBS job re-runable
mv ${jname}.key ${jname}.key.done


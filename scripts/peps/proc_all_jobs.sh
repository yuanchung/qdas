#!/bin/sh
#
# Processing all un-packed rc peps jobs
#

# modify the following line for correct job names
PREFIX="rcox"
LIST=`find . -name rcox_000_0000.out`

# no need to change after this line
basedir=`pwd`

for i in $LIST
do
  dname=`dirname ${i}`
  echo "Checking ${dname}..."
  cd ${dname}
  nkey=`ls ${PREFIX}*.key ${PREFIX}*.key.done | wc -l`
  npbslog=`grep 'Int. Signal' ${PREFIX}*.out | wc -l`
  # check whether the job is finished
  if [ "$nkey" = "$npbslog" ] ; then
    # a complete job, process it!
    echo " complete job!"
    echo " Processing ${dname}..."
    # change all .key.done to .key
    for j in ${PREFIX}*.key.done
    do
      bname=`basename ${j} .done`
      mv ${j} ${bname}
    done
    proc_3ppeps.sh ${PREFIX} > /dev/null
    tar cfz data.tgz ${PREFIX}*.key ${PREFIX}*.out
    rm ${PREFIX}*.key ${PREFIX}*.out
  else
    echo "  not completed; leave it along."
  fi
  cd $basedir
done


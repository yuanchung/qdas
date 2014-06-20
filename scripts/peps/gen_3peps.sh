#!/bin/sh
# generate a bunch of 3PPES jobs

# the prefix of template files; two files must be presented
# in the working directory:
#  ${JOBNAME}.key-tmpl -> the template keyword file
#  ${JOBNAME}.pbs-tmpl -> the template PBS jobfile
JOBNAME="rcox"

# base time for the first pulse
TAU0=100
TLIST="0 10 20 30 40 50 60 70 80 90 100 125 150 175 200 300 400 500 750 1000"
TAU_MIN=-12
TAU_MAX=36
TAU_STEP=3

# generate pbs file
dir=`pwd | sed 's/\/.*\///g'`
echo "Generating PBS script: ${dir}.pbs..."
cat ${JOBNAME}.pbs-tmpl | sed "s/JOBNAME/${dir}/g" > ${dir}.pbs

for T in $TLIST
do
  T_str=`printf "%04d" ${T}`
  tau=$TAU_MIN
  while [ $tau -le $TAU_MAX ] 
  do
    tau_str=`printf "%03d" ${tau}`
    t1=0
    t2=$tau
    t3=`expr $tau + ${T}`
    # shift {t1,t2,t3} so that the smallest time is TAU0
    smallest=`printf "%d\n%d\n%d\n" ${t1} ${t2} ${t3} | sort -n | head -n 1`
    t1=`expr ${t1} + ${TAU0} - ${smallest}`
    t2=`expr ${t2} + ${TAU0} - ${smallest}`
    t3=`expr ${t3} + ${TAU0} - ${smallest}`
    echo "Generating ${JOBNAME}_${tau_str}_${T_str}.key"
    cat ${JOBNAME}.key-tmpl | sed "s/TAU1/${t1}/g" |\
            sed "s/TAU2/${t2}/g" |\
            sed "s/TAU3/${t3}/g" > ${JOBNAME}_${tau_str}_${T_str}.key

    tau=`expr $tau + $TAU_STEP`
  done
done

# DONE


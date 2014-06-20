#!/bin/sh
# generate a bunch of 3PPES jobs for 2D spectrum

# with positive and negative tau
# note that the negative tau jobs used type-2 scane sequence,
# i.e. non-rephasing sequence corresponding to signal at the
# k1-k2+k3 direction

# the prefix of template files; two files must be presented
# in the working directory:
#  ${JOBNAME}.key-tmpl -> the template keyword file
#  ${JOBNAME}.pbs-tmpl -> the template PBS jobfile
JOBNAME="three"

# base time for the first pulse
TAU0=80
TLIST="0 100 300"
# TAU_MAX must be positive
TAU_MAX=300
TAU_STEP=10

dir=`pwd | sed 's/\/.*\///g'`

for T in $TLIST
do
  T_str=`printf "%04d" ${T}`
  tau=0
  while [ $tau -le $TAU_MAX ] 
  do
    # positive tau
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

    # negative tau; pulse2 ahead of pulse1
    if [ "${tau}" != "0" ] ; then
      t1=$tau
      t2=0
      # T is still defined as the time period between 
      # the second and third pulses
      t3=`expr $tau + ${T}`
      # shift {t1,t2,t3} so that the smallest time is TAU0
      smallest=`printf "%d\n%d\n%d\n" ${t1} ${t2} ${t3} | sort -n | head -n 1`
      t1=`expr ${t1} + ${TAU0} - ${smallest}`
      t2=`expr ${t2} + ${TAU0} - ${smallest}`
      t3=`expr ${t3} + ${TAU0} - ${smallest}`
      echo "Generating ${JOBNAME}_-${tau_str}_${T_str}.key"
      cat ${JOBNAME}.key-tmpl | sed "s/TAU1/${t1}/g" |\
              sed "s/TAU2/${t2}/g" |\
              sed "s/TAU3/${t3}/g" > ${JOBNAME}_-${tau_str}_${T_str}.key
    fi
    tau=`expr $tau + $TAU_STEP`
    sleep 0.5
  done
done

# DONE


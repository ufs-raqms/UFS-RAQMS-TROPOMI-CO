#!/bin/bash -x
export REMAPSHTR=$FV3DIR_RELEASE/ush/fv3gfs.5degll.only.sh    
export new_date=/home/lenzen/UTILS/da_advance_time.exe
export INC=${INC:-6}
echo INC $inc FHMAX $FHMAX
echo SCRATCH_OUT $SCRATCH_OUT
sleep 5
#for (( i=0 ; i <= $FHMAX ; i += $INC)); do
echo CDATE $CDATE
for (( i=$INC ; i <= $FHMAX ; i += $INC)); do
  export CDATEF=$($new_date $CDATE ${i}h)
  echo CDATEF $CDATEF
  $REMAPSHTR >& $SCRATCH_OUT/outregrid.tracer.$CDATEF   &
done
wait

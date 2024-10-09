#/bin/bash -x
export FIXDYNEXE=$FV3DIR_RELEASE/exec/fixdyn
export REMAPSHTR=$FV3DIR_RELEASE/ush/fv3gfs_tracer.CONTROL.sh    
export REMAPSHDYN=$FV3DIR_RELEASE/ush/fv3gfs_dynam.sh    
#export CDATEF=2019072118
export new_date=/home/lenzen/UTILS/da_advance_time.exe
#echo FHMAx $FHMAX
export INC=6
#for (( i=0 ; i <= $FHMAX ; i += $INC)); do
for (( i=$INC ; i <= $FHMAX ; i += $INC)); do
  export CDATEF=$($new_date $CDATE ${i}h)
#  echo $i $CDATEF
  $REMAPSHTR >& $DATAOUTREGRID/outregrid.tracer.DIAG.$CDATEF  &
  export HHH=`printf "%03d " "$i"`
#export REMAP_LAUNCHER="srun -prepend-rank -n 1"
#  echo $i $HHH
  $REMAPSHDYN >& $DATAOUTREGRID/outdyn.CONTROL.$CDATEF &
done
wait

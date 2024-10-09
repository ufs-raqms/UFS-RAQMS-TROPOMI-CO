#!/bin/bash -x
export CDATE=2019093000
export CDATE=2019123118
#export CDATEEND=2019100112
export CDATEEND=2019100106
export CDATEEND=2020010106
export SCRATCH_OUT=/scratch/users/lenzen/OUTTR
mkdir $SCRATCH_OUT
export FHMAX=6
export FV3DIR_RELEASE=/home/lenzen/EMC_FV3/v91/EMC_FV3GFS-GSDCHEM/release/v0/
#export BASE_PATH=/ships19/models2/lenzen/FV3GFS.9.1.2019/O3.CONTROL.CEDS.run/C192/
export BASE_PATH=/ships22/raqms/lenzen/FV3GFS.9.1.2019/O3.CONTROL.CEDS.CONT.run/C192/
export BASE_OUT=$BASE_PATH
export LL_OUT=/ships22/raqms/lenzen/FV3GFS.9.1.2019/O3.CONTROL.CEDS.CONT.run/C192/

export REMAPSHTR=$FV3DIR_RELEASE/ush/fv3gfs.5degll.only.sh    
export new_date=/home/lenzen/UTILS/da_advance_time.exe
export INC=${INC:-6}
echo INC $inc FHMAX $FHMAX
echo SCRATCH_OUT $SCRATCH_OUT
sleep 5
#for (( i=0 ; i <= $FHMAX ; i += $INC)); do
echo CDATE $CDATE
while (( CDATE <= CDATEEND )) ;
 do
for (( i=$INC ; i <= $FHMAX ; i += $INC)); do
  export DATATRACE=$BASE_OUT/TRACER/$CDATE/
  export CDATEF=$($new_date $CDATE ${i}h)
  export PATHLL=$LL_OUT/5DEGll/
  echo CDATEF $CDATEF
  $REMAPSHTR >& $SCRATCH_OUT/outregrid.tracer.$CDATEF   &
  export CDATE=$CDATEF
  sleep 5
done
done
wait

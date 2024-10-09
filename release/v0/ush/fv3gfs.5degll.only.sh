#!/bin/bash -x

export CDATE=${CDATE:-"2016100300"}
export CASE=${CASE:-"C192"}           # C48 C96 C192 C384 C768 C1152 C3072

pwd=$(pwd)
export DATARE=${DATATRACE:-$pwd}

cd $DATARE || exit 8

export FULLTRACER=FULL
# DATAREGRID in output directory from above
# PATHP5DEG is output directory for below
echo DATAREGRID $DATAREGRID
echo at $pwd
pwd
export DATAFV3IN=$DATATRACE
export PATHLL=${PATHLL:-$BASE_PATH/5DEGll/}
mkdir -p $PATHLL
echo do fv32ll.bilin.uvpole.x CDATEF $CDATEF
srun -n 1 /home/lenzen/UTILS/fv32ll.bilin.uvpole.x >$SCRATCH_OUT/outfv32ll.bilin.uvpole.$CDATEF

exit $err

#!/bin/bash -x
/bin/rm -rf $PATHSCRATCH/input.ready
export PPIDOLD=${PPID}
export GSI_DIR=$GSI_DIR/sorc/gsi.fd/
sleep 120
export CDATE0=$CDATE
mkdir -p $PATHSCRATCH
export DONE=false
export gsit=0
until [[ "$DONE" == true ]]
  do
gsit=$((gsit+1))
times=0
echo PATHSCRATCH $PATHSCRATCH
maxtime=3600
until [[ -s $PATHSCRATCH/input.ready || $times -gt $maxtime ]]
  do
    sleep 30
    times=$((times+30))
    echo times $times
    export PPIDNEW=$PPID
    if [[ ! $PPIDNEW == $PPIDOLD ]] ; then
        echo parent dead 
        exit -1
    fi  
    echo PPID $PPID
  done
echo seconds wait $times gsit $gsit CDATE $CDATE
source $PATHSCRATCH/cdate.inc
export gsi_scratch=$SCRATCH_OUT/GSISCRATCH/$CDATE0
export gsi_path=$BASE_OUT/GSIDIR/$CDATE0
export GSIDIR=$gsi_path
mkdir -p $gsi_scratch
export gsi_scratch0=$gsi_scratch
cp -p $PATHSCRATCH/cdate.inc $PATHSCRATCH/cdate.inc.laSt
cp -p $PATHSCRATCH/input.ready $PATHSCRATCH/input.ready.last
/bin/rm -f $PATHSCRATCH/input.ready $PATHSCRATCH/cdate.inc
echo CDATE $CDATE
#Eexport DATAREGRID=$PWD/../../GSIDIR/$CDATE0/
export DATAREGRID=$BASE_PATH/GSIDIR/$CDATE0/
#export DATAINC=$PWD/../../GSIINC/$CDATE0/
export DATAINC=$BASE_PATH/GSIINC/$CDATE0/
mkdir -p $DATAINC
export CASE=C192
#mkdir -p $DATAREGRID
export BASE_DATA=/ships19/aqda/lenzen/FV3GFS_V1_RELEASE 
#export REMAPDIR=/home/lenzen/FV3GFS-GSDCHEM.0.9.0/EMC_FV3GFS-GSDCHEM/release/v0/
#export REMAPEXE=${REMAPEXE:-$REMAPDIR/exec/fregrid_parallel}
export REMAP_LAUNCHER="srun -s -n 11 "
export FIX_FV3=$BASE_DATA/fix/fix_fv3  
export grid_loc=$FIX_FV3/$CASE/${CASE}_mosaic.nc
#export DATATRACE=/ships19/aqda/lenzen/RAQMSEMISNC/C192/
export DATA=$DATAREGRID
export master_grid="0p5deg"
export GG=${master_grid:-"0p25deg"}   # 1deg 0p5deg 0p25deg 0p125deg
export weight_file=$FIX_FV3/$CASE/remap_weights_${CASE}_${GG}.nc
export tracer="ps, zsurf, sphum, delp, pres, T, U, V, co, incco"
export CDATEF=$CDATE
export DATAFV3IN=$DATAREGRID
export PATHLL=$DATAREGRID
sleep 2 # to stop failure
#/home/lenzen/UTILS/fv32ll.x >& /scratch/users/lenzen/GSISCRATCH/out.$CDATE.outfv32ll
echo SLURM_JOB_NUM_NODE $SLURM_JOB_NUM_NODE
srun -O --mem=6G  -n 1 /home/lenzen/UTILS/fv32ll.gen.gen.deflate.bilin.uvpole.x >& $gsi_scratch/out.$CDATE.outfv32ll.deflate
if [ $GG = 1deg    ];  then  export nlon=360  ; export nlat=180  ; fi
if [ $GG = 0p5deg  ];  then  export nlon=720  ; export nlat=360  ; fi
if [ $GG = 0p25deg ];  then  export nlon=1440 ; export nlat=720  ; fi
if [ $GG = 0p125deg ]; then  export nlon=2880 ; export nlat=1440 ; fi
err=0
export CDATEF=$CDATE
#export DATASCRATCH=/scratch/users/lenzen/GSISCRATCH/$EXPNAME/$CDATE0/$CDATEF/
export DATASCRATCH=$gsi_scratch/$EXPNAME/$CDATE0/$CDATEF/
echo DATASCRATCH $DATASCRATCH
mkdir -p $DATASCRATCH
export scratch_d_gsi=$DATASCRATCH
export GSIDATE=$CDATEF
export BUFRDATE=$CDATEF
export GSIDIR=$gsi_path
#export GSI_DIR=/home/lenzen/GSI/FV3GFS.TROPOMI/
env | grep NCPUS
echo finished reamap >& $gsi_scratch/outremapdone.$CDATEF
echo GSI_DIR $GSI_DIR/run/rungsi.FV3gfs.GSI.CO.sh >& $gsi_scratch/outgsirun.$CDATEF
$GSI_DIR/run/run_gsi_TROPOMI.CO.ak.blended.cofire.4xtime.ksh >& $gsi_scratch/outgsi.CO.TROPOMI.blended.$CDATEF
export PATH_GSIINC=$DATAINC
export gsi_scratch=$scratch_d_gsi
rsync -ltvD $gsi_scratch/sigf00 $PATH_GSIINC/sigf00.$CDATEF
/home/lenzen/UTILS/maplltofv3.x >&$gsi_scratch/outmaplltofv3.$CDATEF
#echo $CDATE >$gsi_scratch/../../gsi.done.$CDATE
echo $CDATE >$PATHSCRATCH/gsi.done.$CDATE
echo gsit $gsit CDATE $CDATE
# ajl new
#export CDATE0=$CDATEF
#export PATHSCRATCH=/scratch/users/lenzen/GSISCRATCH/$EXPNAME/$CDATE0/
#if [ $gsit -ge 9 ] ; then
#$  export DONE=true
#$fi
done

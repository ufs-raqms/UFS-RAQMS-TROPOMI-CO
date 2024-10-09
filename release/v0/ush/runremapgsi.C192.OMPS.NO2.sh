#!/bin/bash -x
/bin/rm $PATHSCRATCH/input.ready
sleep 120
#env >/scratch/users/lenzen/outenv.gsi
echo gsi_path $gsi_path
echo path_bsi $path_gsi
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
until [[ -s $PATHSCRATCH/input.ready || $times -gt $maxtime  ]]
  do
    sleep 30
    times=$((times+30))
    echo times $times
  done
echo seconds wait $times gsit $gsit CDATE $CDATE
source $PATHSCRATCH/cdate.inc
cat $PATHSCRATCH/cdate.inc
cp -p $PATHSCRATCH/cdate.inc $PATHSCRATCH/cdate.inc.laSt
cp -p $PATHSCRATCH/input.ready $PATHSCRATCH/input.ready.last
/bin/rm -f $PATHSCRATCH/input.ready $PATHSCRATCH/cdate.inc
echo CDATE $CDATE
echo CDATE $CDATE times $times >$gsi_scratch/outcdate.$CDATE
export DATAREGRID=$BASE_PATH/GSIDIR/$CDATE0/
mkdir -p $DATAREGRID
export DATAINC=$BASE_PATH/GSIINC/$CDATE0/
mkdir -p $DATAINC
export CASE=C192
export BASE_DATA=/ships19/aqda/lenzen/FV3GFS_V1_RELEASE 
export REMAPDIR=/home/lenzen/FV3GFS-GSDCHEM.0.9.0/EMC_FV3GFS-GSDCHEM/release/v0/
export REMAPEXE=${REMAPEXE:-$REMAPDIR/exec/fregrid_parallel}
export REMAP_LAUNCHER="srun -s --mem=10G -n 11 "
export FIX_FV3=$BASE_DATA/fix/fix_fv3  
export grid_loc=$FIX_FV3/$CASE/${CASE}_mosaic.nc
export DATA=$DATAREGRID
export master_grid="0p5deg"
export GG=${master_grid:-"0p25deg"}   # 1deg 0p5deg 0p25deg 0p125deg
export weight_file=$FIX_FV3/$CASE/remap_weights_${CASE}_${GG}.nc
export tracer="ps, zsurf, pblht, ptrop, sphum, delp, pres, T, U, V, no2,incno2"
export CDATEF=$CDATE
export DATAFV3IN=$DATAREGRID
export PATHLL=$DATAREGRID
#srun -O --mem=6G  -n 1 /home/lenzen/UTILS/fv32ll.gen.gen.deflate.x >& $gsi_scratch/out.$CDATE.outfv32ll.deflate.$CDATE
sleep 3
srun -O --mem=6G  -n 1 /home/lenzen/UTILS/fv32ll.bilin.uvpole.x >& $gsi_scratch/out.$CDATE.outfv32ll.deflate.$CDATE
sleep 3
export tracer="ps, zsurf, sphum, delp, pres, T, U, V, no2,incno2"
if [ $GG = 1deg    ];  then  export nlon=360  ; export nlat=180  ; fi
if [ $GG = 0p5deg  ];  then  export nlon=720  ; export nlat=360  ; fi
if [ $GG = 0p25deg ];  then  export nlon=1440 ; export nlat=720  ; fi
if [ $GG = 0p125deg ]; then  export nlon=2880 ; export nlat=1440 ; fi
err=0
export DATASCRATCH=$SCRATCH_OUT/$CDATE0/$CDATEF/
echo DATASCRATCH $DATASCRATCH
mkdir -p $DATASCRATCH
export scratch_d_gsi=$DATASCRATCH
export GSIDATE=$CDATEF
export BUFRDATE=$CDATEF
env | grep NCPUS
echo finished reamap >& $SCRATCH_OUT/outremapdone.$CDATEF
echo GSI_DIR $GSI_DIR/run/run_gsi_FV3GFS.omps.no2.ksh  >& $SCRATCH_OUT/outgsirun.omps.no2.$CDATEF
#env > /scratch/users/lenzen/out.before.rungsi
#echo gsi_path before $gsi_path
export OMPSNO2PATH=/ships19/aqda/lenzen/OMPSNO2NC/
#env >/scratch/users/lenzen/outenv.before.run.gsi2
echo GSIPROC $GSIPROC
echo GSIPROCNO2 $GSIPROCNO2
echo SLURM_NTASKS $SLURM_NTASKS
echo NTASKS $NTASKS
echo SLURM_NTASKS $SLURM_NTASKS
$GSI_DIR/sorc/gsi.fd/run/run_gsi_FV3GFS.omps.no2.ksh >& $SCRATCH_OUT/outrun_gsi.OMPS.NO2.TROP.$CDATEF
export PATH_GSIINC=$DATAINC
export gsi_scratch=$scratch_d_gsi
sleep 3
rsync -ltvD $gsi_scratch/sigf00 $PATH_GSIINC/sigf00.$CDATEF
#rsync -ltvD $gsi_scratch/ges.nc $PATH_GSIINC/ges.$CDATEF.nc
/home/lenzen/UTILS/maplltofv3.x >&$gsi_scratch/outmaplltofv3.$CDATEF
sleep 3
echo $CDATE >$PATHSCRATCH/gsi.done.$CDATE
echo gsit $gsit CDATE $CDATE
done

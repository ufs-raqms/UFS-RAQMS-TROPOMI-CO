#!/bin/bash  -x
export CASE=C192
export CDATE=2019071512
export GG="0p5deg"
export BASE_DATA=/home/lenzen/FV3GFS-GSDCHEM.0.9.0.VIIRS.AOD/EMC_FV3GFS-GSDCHEM/release/v0/
export REMAPEXE=${REMAPEXE:-$BASE_DATA/exec/fregrid_parallel}
export BASE_DATA=/ships19/aqda/lenzen/FV3GFS_V1_RELEASE ;# data directory
export INDIR=/ships19/models2/lenzen/FV3GFS.9.1.2019/O3.CONTROL.CEDS.run/C192/TRACER/2019071512/
export FIX_FV3=${FIX_FV3:-$BASE_DATA/fix/fix_fv3}
export grid_loc=$FIX_FV3/$CASE/${CASE}_mosaic.nc
export weight_file=$FIX_FV3/$CASE/remap_weights_${CASE}_${GG}.nc
export APRRUN="srun -n 1"

export REMAP_LAUNCHER=${REMAP_LAUNCHER:-${APRUN:-""}}
export NTHREADS_REMAP=${NTHREADS_REMAP:-${NTHREADS:-1}}

#--------------------------------------------------
if [ $GG = 1deg    ];  then  export nlon=360  ; export nlat=180  ; fi
if [ $GG = 0p5deg  ];  then  export nlon=720  ; export nlat=360  ; fi
if [ $GG = 0p25deg ];  then  export nlon=1440 ; export nlat=720  ; fi
if [ $GG = 0p125deg ]; then  export nlon=2880 ; export nlat=1440 ; fi


export INDIR=/ships19/models2/lenzen/FV3GFS.9.1.2019/O3.CONTROL.CEDS.run/C192/TRACER/2019071512/
export DATA=$INDIR
export DATAREGRID=/ships22/raqms/lenzen/CDES.EMIS
export CDATEF=2019071512
mkdir -p $DATAREGRID
#export srcpropane=srcpropane.nc
export srcethane=srcethane.nc
export srcpropene=srcpropene.nc
export srcpentane=srcpentane.nc
export srcbutane=srcbutane.nc
export srcco_totl=srcco_totl.nc
export srcnind=srcnind.nc
export srchexane=srchexane.nc
export types="srcbutane"
for type in $types ; do
  echo type $type
  export in_file="${type}"

  export fld=$(eval echo \${${type}})
  echo in_file $in_file
  echo fld $fld
  export out_file="${type}.${CDATEF}.nc"
  [[ -s $DATAREGRID/$out_file ]] && rm -f $DATAREGRID/$out_file
  export fld=$(eval echo \${${type}})
  echdo in_file $in_file

  $REMAP_LAUNCHER $REMAPEXE --input_dir $DATA \
                         --input_file $in_file \
                         --output_dir $DATAREGRID \
                         --output_file $out_file \
                         --input_mosaic $grid_loc \
                         --scalar_field "$fld" \
                         --interp_method conserve_order1 \
                         --remap_file $weight_file \
                         --nlon $nlon \
                         --debug \
                         --nlat $nlat
  rc=$?
  ((err+=$rc))

done

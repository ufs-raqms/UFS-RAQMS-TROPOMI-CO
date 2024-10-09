#!/bin/bash -x
################################################################################
# UNIX Script Documentation Block
# Script name:         exglobal_fcst_fv3gfs.sh.ecf
# Script description:  Runs a global FV3GFS model forecast
#
# Author:   Fanglin Yang       Org: NCEP/EMC       Date: 2016-11-15
# Abstract: This script runs a single GFS forecast with FV3 dynamical core.
#           This script is created based on a C-shell script that GFDL wrote
#           for the NGGPS Phase-II Dycore Comparison Project.
#
# Script history log:
# 2016-11-15  Fanglin Yang   First Version.
# 2017-02-09  Rahul Mahajan  Added warm start and restructured the code.
# 2017-03-10  Fanglin Yang   Updated for running forecast on Cray.
# 2017-03-24  Fanglin Yang   Updated to use NEMS FV3GFS with IPD4 
#
# Attributes:
#   Language: Portable Operating System Interface (POSIX) Shell
#   Machine: WCOSS-CRAY, Theia
################################################################################

#  Set environment.
ulimit -a >$HOME/outulimit.submit
export VERBOSE="NO"
export VERBOSE=${VERBOSE:-"YES"}
if [ $VERBOSE = YES ] ; then
  echo $(date) EXECUTING $0 $* >&2
  set -x
fi

# Cycling and forecast hour specific parameters
export PSLOT=${PSLOT:-fv3gfs}
export CASE=${CASE:-C768}
export CDATE=${CDATE:-2017032500}
export CDUMP=${CDUMP:-gfs}
export FHMIN=${FHMIN:-0}
export FHMAX=${FHMAX:-240}
export FHOUT=${FHOUT:-3}
export FHZER=${FHZER:-6}
export FHCYC=${FHCYC:-24}

# Directories.
export PTMP=${PTMP:-/gpfs/hps/ptmp}
export STMP=${STMP:-/gpfs/hps/stmp}
export NWPROD=${NWPROD:-${NWROOT:-/nwprod}}
export BASE_DATA=${BASE_DATA:-$NWPROD}
export FIX_DIR=${FIX_DIR:-$BASE_DATA/fix}
export FIX_AM=${FIX_AM:-$FIX_DIR/fix_am}
export FIX_FV3=${FIX_FV3:-$FIX_DIR/fix_fv3}
export DATA=${DATA:-$STMP/$LOGNAME/pr${PSLOT}${CASE}_${CDATE}}    #temporary running directory
export ROTDIR=${ROTDIR:-$PTMP/$LOGNAME/pr${PSLOT}}                #rorating archive directory
export IC_DIR=${IC_DIR:-$PTMP/$LOGNAME/ICs}                       #cold start initial conditions

# Model resolution specific parameters
export DELTIM=${DELTIM:-225}
export layout_x=${layout_x:-8}
export layout_y=${layout_y:-16}
export LEVS=${LEVS:-64}

# Utilities
#export NCP=${NCP:-"/bin/cp -p"}
export NCP="rsync -ltvD"
export NCP=${NCP:-"/bin/cp -p"}
export NLN=${NLN:-"/bin/ln -sf"}
export SEND=${SEND:-"YES"}   #move final result to rotating directory
export ERRSCRIPT=${ERRSCRIPT:-'eval [[ $err = 0 ]]'}
export NDATE=${NDATE:-$NWPROD/util/exec/ndate}

# Other options
export MEMBER=${MEMBER:-"-1"} # -1: control, 0: ensemble mean, >0: ensemble member $MEMBER
export ENS_NUM=${ENS_NUM:-1}

# Model specific stuff
export FCSTEXECDIR=${FCSTEXECDIR:-${EXECDIR:-$BASE_DATA/sorc/fv3gfs.fd/BUILD/bin}}
export FCSTEXEC=${FCSTEXEC:-fv3_gfs.x}
export PARM_FV3DIAG=${PARM_FV3DIAG:-$FV3DIR_RELREASE/parm/parm_fv3diag}

# Model config options
export FCST_LAUNCHER=${FCST_LAUNCHER:-${APRUN:-""}}
export tasks=${tasks:-$((6*layout_x*layout_y))}
export nthreads=${nthreads:-${nth_f:-1}}
export cores_per_node=${cores_per_node:-${task_per_node:-24}}
export ntiles=${ntiles:-6}
export TYPE=${TYPE:-nh}                  # choices:  nh, hydro
export MONO=${MONO:-non-mono}            # choices:  mono, non-mono
export use_hyper_thread=${hyperthread:-".false."}

#-------------------------------------------------------
if [ ! -d $ROTDIR ]; then mkdir -p $ROTDIR; fi
if [ ! -d $DATA ]; then mkdir -p $DATA ;fi
mkdir -p $DATA/RESTART $DATA/INPUT $DATA/GBBEPx
rsync -ltvD $GBBEPxDIR/*C192* $DATA/GBBEPx
mkdir -p $DATA/GBBEPx/Inputfire/
rsync -ltvD /ships19/aqda/lenzen/RAQMSEMIS/C192/Inputfire/* $DATA/GBBEPx/Inputfire/
cd $DATA/GBBEPx
for i in {1..6}
  do
    export ctile=tile$i
    mkdir $ctile
    if [ -f $DATA/GBBEPx/GBBEPx.bc.$YYYYMMDD.FV3.$CASE"Grid."$ctile.bin ] ; then
      $NLN $DATA/GBBEPx/GBBEPx.bc.$YYYYMMDD.FV3.$CASE"Grid."$ctile.bin $ctile/ebu_bc.dat
    else
      $NLN $DATA/GBBEPx/GBBEPx.bc.$YYYYMMDD.FV3.$CASE"Grid."$ctile.nc $ctile/ebu_bc.nc
    fi
    if [ -f $DATA/GBBEPx/GBBEPx.bc.$YYYYMMDD.FV3.$CASE"Grid."$ctile.bin ] ; then
      $NLN $DATA/GBBEPx/GBBEPx.oc.$YYYYMMDD.FV3.$CASE"Grid."$ctile.bin $ctile/ebu_oc.dat
    else
      $NLN $DATA/GBBEPx/GBBEPx.oc.$YYYYMMDD.FV3.$CASE"Grid."$ctile.nc $ctile/ebu_oc.nc
    fi
    if [ -f $DATA/GBBEPx/GBBEPx.so2.$YYYYMMDD.FV3.$CASE"Grid."$ctile.bin ] ; then
      $NLN $DATA/GBBEPx/GBBEPx.so2.$YYYYMMDD.FV3.$CASE"Grid."$ctile.bin $ctile/ebu_so2.dat
    else
      $NLN $DATA/GBBEPx/GBBEPx.so2.$YYYYMMDD.FV3.$CASE"Grid."$ctile.nc $ctile/ebu_so2.nc
    fi
    if [ -f $DATA/GBBEPx/GBBEPx.pm25.$YYYYMMDD.FV3.$CASE"Grid."$ctile.bin ] ; then
      $NLN $DATA/GBBEPx/GBBEPx.pm25.$YYYYMMDD.FV3.$CASE"Grid."$ctile.bin $ctile/ebu_pm_25.dat
    else
      $NLN $DATA/GBBEPx/GBBEPx.pm25.$YYYYMMDD.FV3.$CASE"Grid."$ctile.nc $ctile/ebu_pm_25.nc
    fi
    if [ -f $DATA/GBBEPx/meanFRP.$YYYYMMDD.FV3.$CASE"Grid."$ctile.bin ] ; then
      $NLN $DATA/GBBEPx/meanFRP.$YYYYMMDD.FV3.$CASE"Grid."$ctile.bin $ctile/plumefrp.dat
    else
      $NLN $DATA/GBBEPx/meanFRP.$YYYYMMDD.FV3.$CASE"Grid."$ctile.nc $ctile/plumefrp.nc
    fi
    $NLN $DATA/GBBEPx/meanFRP.$YYYYMMDD.FV3.$CASE"Grid."$ctile.nc plumefrp.$ctile.nc
    $NLN $DATA/GBBEPx/meanFRP.$YYYYMMDD.FV3.$CASE"Grid."$ctile.nc $ctile/plumefrp.nc
    $NLN $DATA/GBBEPx/GBBEPx.co.$YYYYMMDD.FV3.$CASE"Grid."$ctile.nc $ctile/ebu_co.nc
    $NLN $DATA/GBBEPx/GBBEPx.co.$YYYYMMDD.FV3.$CASE"Grid."$ctile.nc ebu_co.$ctile.nc
    
  done
cd $DATA || exit 8

#-------------------------------------------------------
# member directory
if [ $MEMBER -lt 0 ]; then
  PREINP=$CDUMP
  MEMCHAR=""
else
  PREINP=enkf.$CDUMP
  MEMCHAR=mem`printf %03i $MEMBER`
fi
yyyymmdd=`echo $CDATE | cut -c1-8`
hh=`echo       $CDATE | cut -c9-10`
MEMDIR=$ROTDIR/${PREINP}.$yyyymmdd/$hh/$MEMCHAR
if [ ! -d $MEMDIR ]; then mkdir -p $MEMDIR; fi

#-------------------------------------------------------
# initial conditions
export warm_start=${warm_start:-".false."}
if [ $warm_start = ".false." ]; then
  $NCP $BASE_DATA/ICs/C192_2018081400/gfs_ctrl.nc $DATA/INPUT/.
 if [ -d $IC_DIR/${CASE}_$CDATE ]; then
  $NCP $IC_DIR/${CASE}_$CDATE/* $DATA/INPUT/.
 else
  for file in $MEMDIR/INPUT/*.nc; do
    file2=$(echo $(basename $file))
    fsuf=`echo $file2 | cut -c1-3`
    if [ $fsuf = "gfs" -o $fsuf = "sfc" ]; then
      $NLN $file $DATA/INPUT/$file2
    fi
  done
 fi
else
  for file in $MEMDIR/RESTART/*.nc; do
    file2=$(echo $(basename $file))
    $NLN $file $DATA/IrNPUT/$file2
  done
  $NLN $MEMDIR/RESTART/coupler.res $DATA/INPUT/coupler.res
  export read_increment=${read_increment:-".false."}
  if [ $read_increment == ".true." ]; then
    if [ -f $MEMDIR/$increment_file ]; then
      $NLN $MEMDIR/$increment_file $DATA/INPUT/$increment_file
    else
      export read_increment=".false."
    fi
  fi
fi
nfiles=`ls -1 $DATA/INPUT/* | wc -l`
if [ $nfiles -lt 0 ]; then
  echo "Initial conditions must exist in $DATA/INPUT, ABORT!"
  exit 1
fi

#--------------------------------------------------------------------------
# Grid and orography data
for n in `seq 1 $ntiles`; do
  $NLN $FIX_FV3/$CASE/${CASE}_grid.tile${n}.nc     $DATA/INPUT/${CASE}_grid.tile${n}.nc
  $NLN $FIX_FV3/$CASE/${CASE}_oro_data.tile${n}.nc $DATA/INPUT/oro_data.tile${n}.nc
  $NLN $FIX_FV3/$CASE/angle.tile${n}.nc $DATA/INPUT/angle.tile${n}.nc
done
$NLN $FIX_FV3/$CASE/${CASE}_mosaic.nc  $DATA/INPUT/grid_spec.nc
# new ajl
$NLN $FIX_FV3/$CASE/${CASE}_mosaic.nc  $DATA/INPUT/atmos_mosaic.nc

# GFS standard input data
export iems=${iems:-1}
export isol=${isol:-2}
export iaer=${iaer:-111}
export ico2=${ico2:-2}

$NLN $FIX_AM/global_solarconstant_noaa_an.txt  $DATA/solarconstant_noaa_an.txt
#$NLN $FIX_AM/global_o3prdlos.f77               $DATA/INPUT/global_o3prdlos.f77
$NLN $FIX_AM/global_o3prdlos.f77               $DATA/global_o3prdlos.f77
$NLN $FIX_AM/global_sfc_emissivity_idx.txt     $DATA/sfc_emissivity_idx.txt

$NLN $FIX_AM/global_co2historicaldata_glob.txt $DATA/co2historicaldata_glob.txt
$NLN $FIX_AM/co2monthlycyc.txt                 $DATA/co2monthlycyc.txt
if [ $ico2 -gt 0 ]; then
  for file in `ls $FIX_AM/fix_co2_proj/global_co2historicaldata* ` ; do
    $NLN $file $DATA/$(echo $(basename $file) | sed -e "s/global_//g")
  done
fi

$NLN $FIX_AM/global_climaeropac_global.txt     $DATA/aerosol.dat
if [ $iaer -gt 0 ] ; then
  for file in `ls $FIX_AM/global_volcanic_aerosols* ` ; do
    $NLN $file $DATA/$(echo $(basename $file) | sed -e "s/global_//g")
  done
fi

export FNGLAC=${FNGLAC:-"$FIX_AM/global_glacier.2x2.grb"}
export FNMXIC=${FNMXIC:-"$FIX_AM/global_maxice.2x2.grb"}
export FNTSFC=${FNTSFC:-"$FIX_AM/RTGSST.1982.2012.monthly.clim.grb"}
export FNSNOC=${FNSNOC:-"$FIX_AM/global_snoclim.1.875.grb"}
export FNZORC=${FNZORC:-"igbp"}
export FNALBC=${FNALBC:-"$FIX_AM/global_snowfree_albedo.bosu.t1534.3072.1536.rg.grb"}
export FNALBC2=${FNALBC2:-"$FIX_AM/global_albedo4.1x1.grb"}
export FNAISC=${FNAISC:-"$FIX_AM/CFSR.SEAICE.1982.2012.monthly.clim.grb"}
export FNTG3C=${FNTG3C:-"$FIX_AM/global_tg3clim.2.6x1.5.grb"}
export FNVEGC=${FNVEGC:-"$FIX_AM/global_vegfrac.0.144.decpercent.grb"}
export FNVETC=${FNVETC:-"$FIX_AM/global_vegtype.igbp.t1534.3072.1536.rg.grb"}
export FNSOTC=${FNSOTC:-"$FIX_AM/global_soiltype.statsgo.t1534.3072.1536.rg.grb"}
export FNSMCC=${FNSMCC:-"$FIX_AM/global_soilmgldas.t1534.3072.1536.grb"}
export FNMSKH=${FNMSKH:-"$FIX_AM/seaice_newland.grb"}
export FNVMNC=${FNVMNC:-"$FIX_AM/global_shdmin.0.144x0.144.grb"}
export FNVMXC=${FNVMXC:-"$FIX_AM/global_shdmax.0.144x0.144.grb"}
export FNSLPC=${FNSLPC:-"$FIX_AM/global_slope.1x1.grb"}
export FNABSC=${FNABSC:-"$FIX_AM/global_mxsnoalb.uariz.t1534.3072.1536.rg.grb"}

# nstf_name contains the NSST related parameters
# nstf_name(1) : 0 = NSSTM off, 1 = NSSTM on but uncoupled, 2 = NSSTM on and coupled
# nstf_name(2) : 1 = NSSTM spin up on, 0 = NSSTM spin up off
# nstf_name(3) : 1 = NSST analysis on, 0 = NSSTM analysis off
# nstf_name(4) : zsea1 in mm
# nstf_name(5) : zsea2 in mm
# nst_anl      : .true. or .false., NSST analysis over lake                       
export nstf_name=${nstf_name:-"2,0,1,0,5"}


#------------------------------------------------------------------
# changeable parameters
# dycore definitions
res=`echo $CASE |cut -c2-5`
resp=`expr $res + 1 `
export npx=$resp
export npy=$resp
export npz=`expr $LEVS - 1 `
export io_layout="1,1"
export ncld=5
#export ncols=$(( (${npx}- 1)*(${npy}-1)*3/2 ))

# blocking factor used for threading and general physics performance
#export nyblocks=`expr \( $npy - 1 \) \/ $layout_y `
#export nxblocks=`expr \( $npx - 1 \) \/ $layout_x \/ 32`
#if [ $nxblocks -le 0 ]; then export nxblocks=1 ; fi
export blocksize=${blocksize:-32}

# export the pre-conditioning of the solution
# =0 implies no pre-conditioning
# >0 means new adiabatic pre-conditioning
# <0 means older adiabatic pre-conditioning
export na_init=${na_init:-1}

# variables for controlling initialization of NCEP/NGGPS ICs
export filtered_terrain=${filtered_terrain:-".true."}
export gfs_dwinds=${gfs_dwinds:-".true."}

# determines whether FV3 or GFS physics calculate geopotential
export gfs_phil=${gfs_phil:-".false."}

# determine whether ozone production occurs in GFS physics
export ozcalc=${ozcalc:-".true."}

# export various debug options
export no_dycore=${no_dycore:-".false."}
export dycore_only=${adiabatic:-".false."}
export chksum_debug=${chksum_debug:-".false."}
# export chksum_debug=".true."
export print_freq=${print_freq:-6}

if [ ${TYPE} = "nh" ]; then
  # non-hydrostatic options
  export make_nh=".true."
  export hydrostatic=".false."
  export phys_hydrostatic=".false."     # can be tested
  export use_hydro_pressure=".false."   # can be tested
  export consv_te="1."
else
  # hydrostatic options
  export make_nh=".false."
  export hydrostatic=".true."
  export phys_hydrostatic=".false."     # will be ignored in hydro mode
  export use_hydro_pressure=".true."    # have to be .true. in hydro mode
  export consv_te="0."
fi

# time step parameters in FV3
export k_split=1
export n_split=8
export n_sponge=30
export cal_pre=.false.
export random_clds=.false.

if [ ${MONO} = "mono" -o ${MONO} = "monotonic" ];  then
  # monotonic options
  export d_con="1."
  export do_vort_damp=".false."
  if [ ${TYPE} = "nh" ]; then
    # non-hydrostatic
    export hord_mt="10"
    export hord_xx="10"
  else
    # hydrostatic
    export hord_mt="10"
    export hord_xx="10"
  fi
else
  # non-monotonic options
  export d_con="1."
  export do_vort_damp=".true."
  if [ ${TYPE} = "nh" ]; then
    # non-hydrostatic
    export hord_mt="6"
    export hord_xx="6"
  else
    # hydrostatic
    export hord_mt="10"
    export hord_xx="10"
  fi
fi

if [ ${MONO} = "non-mono" -a ${TYPE} = "nh" ]; then
  export vtdm4="0.02"
else
  export vtdm4="0.05"
fi

if [ $warm_start = ".false." ]; then # CHGRES'd GFS analyses
  export external_ic=".true."
  export mountain=".false."
  export read_increment=".false."
  export res_latlon_dynamics='""'
else # warm start from restart file
  export external_ic=".false."
  export mountain=".true."
  export make_nh=".false."
  export na_init=0                
  if [ $read_increment = ".true." ]; then # add increments on the fly to the restarts
    export res_latlon_dynamics="$increment_file"
  else
    export res_latlon_dynamics='""'
  fi
fi

# build the date for curr_date and diag_table from CDATE
export SYEAR=`echo $CDATE | cut -c1-4`
export SMONTH=`echo $CDATE | cut -c5-6`
export SDAY=`echo $CDATE | cut -c7-8`
export SHOUR=`echo $CDATE | cut -c9-10`
export YYYYMMDD=`echo $CDATE | cut -c1-8`
export curr_date="${SYEAR},${SMONTH},${SDAY},${SHOUR},0,0"
export restart_secs=${restart_secs:-0}

# copy over the tables
export DIAG_TABLE=${DIAG_TABLE:-$PARM_FV3DIAG/diag_table}
export DATA_TABLE=${DATA_TABLE:-$PARM_FV3DIAG/data_table}
export FIELD_TABLE=${FIELD_TABLE:-$PARM_FV3DIAG/field_table}
echo cplflx $cplflx
#export NHTUW=12
#export NHTUW=2

# build the diag_table with the experiment name and date stamp
cat > diag_table << EOF
FV3 Forecast
$SYEAR $SMONTH $SDAY $SHOUR 0 0
EOF
cat $DIAG_TABLE >> diag_table
$NCP $AODLUT/* .
$NCP $DATA_TABLE data_table
echo FIELD_TABLE $FIELD_TABLE
$NCP $FIELD_TABLE field_table
cat field_table
#ls -l field_table

#pwd >& $HOME/outpwdtable
#ls -l field_table >& $HOME/outlsfieldtable
echo cat nems.configure
#------------------------------------------------------------------
cat > nems.configure <<EOF
#############################################
####  NEMS Run-Time Configuration File  #####
#############################################

# EARTH #
EARTH_component_list: ATM CHM RCHM
EARTH_attributes::
  Verbosity = 0 
::

# ATM #
ATM_model:                      fv3 
ATM_petlist_bounds:             -1 -1
ATM_attributes::
  Verbosity = 0 
::

# CHM #
CHM_model:                      gsdchem
CHM_petlist_bounds:             -1 -1
CHM_attributes::
  Verbosity = 0 
::

# RCHM #
RCHM_model:                      raqmschem
RCHM_petlist_bounds:             -1 -1
RCHM_attributes::
  Verbosity = 0 
::

# Run Sequence #
runSeq::
  @450
    ATM phase1
    ATM -> CHM 
    ATM -> RCHM 
    CHM 
    RCHM 
    CHM -> ATM  
    RCHM -> ATM  
    ATM phase2
  @
::

EOF

export quilting=.true.
#export quilting=.false.
export write_groups=1
echo quilting $quilting
echo cat model_configure
cat > model_configure <<EOF
  total_member:            $ENS_NUM
  PE_MEMBER01:             $tasks
  start_year:              $SYEAR
  start_month:             $SMONTH
  start_day:               $SDAY
  start_hour:              $SHOUR
  start_minute:            0
  start_second:            0
  nhours_fcst:             $FHMAX
  RUN_CONTINUE:            ${RUN_CONTINUE:-".false."}
  ENS_SPS:                 ${ENS_SPS:-".false."}

  dt_atmos:                $DELTIM    
  calendar:                ${calendar:-'julian'}
  memuse_verbose:          ${memuse_verbose:-".false."}
  atmos_nthreads:          $nthreads
  use_hyper_thread:        ${hyperthread:-".false."}
  ncores_per_node:         $cores_per_node
  restart_interval:        ${restart_interval:-0}
  quilting:                ${quilting:-".false."}
  write_groups:            ${write_groups:-0}
  write_tasks_per_group:   6
  num_files:               2
  filename_base:          'dyn' 'phy'
  output_grid:             "cubed_sphere_grid"
  nfhout:                  6
  nfhmax_hf:               12
  nfhout_hf:                3
  nsout:                   -1
  cpl:                     .true.
  atm_coupling_interval_sec:  $DELTIM
EOF
cat model_configure > /home/lenzen/outmodelconfig

#&coupler_nml
#  months = ${months:-0}
#  days = ${days:-$((FHMAX/24))}
#  hours = ${hours:-$((FHMAX-24*(FHMAX/24)))}
#  dt_atmos = $DELTIM
#  dt_ocean = $DELTIM
#  current_date = $curr_date
#  calendar = 'julian'
#  memuse_verbose = .false.
#  atmos_nthreads = $nthreads
#  use_hyper_thread = ${hyperthread:-".false."}
#  ncores_per_node = $cores_per_node
#  restart_secs = $restart_secs   ##DA
#/


echo cat input.nml >/home/lenzen/catinput.nml
cat > input.nml <<EOF
&amip_interp_nml
  interp_oi_sst = .true.
  use_ncep_sst = .true.
  use_ncep_ice = .false.
  no_anom_sst = .false.
  data_set = 'reynolds_oi'
  date_out_of_range = 'climo'
/

&atmos_model_nml
  blocksize = $blocksize
  chksum_debug = $chksum_debug
  dycore_only = $dycore_only
  fdiag = ${fdiag:-$FHOUT}
/

&diag_manager_nml
  prepend_date = .F.   
/

&fms_io_nml
  checksum_required = .false.
  max_files_r = 100
  max_files_w = 100
/

&fms_nml
  clock_grain = 'ROUTINE'
  domains_stack_size = ${domains_stack_size:-115200}
  print_memory_usage = ${print_memory_usage:-".false."}
/

&fv_grid_nml
  grid_file = 'INPUT/grid_spec.nc'
/

&fv_core_nml
  layout = $layout_x,$layout_y
  io_layout = $io_layout
  npx = $npx
  npy = $npy
  ntiles = $ntiles
  npz = $npz
  grid_type = -1
  make_nh = $make_nh
  fv_debug = ${fv_debug:-".false."}
  range_warn = ${range_warn:-".false."}
  reset_eta = .false.
  n_sponge = ${n_sponge:-24}
  nudge_qv = ${nudge_qv:-".true."}
  rf_fast = .false.
  tau = 5.
  rf_cutoff = 7.5e2
  d2_bg_k1 = 0.15
  d2_bg_k2 = 0.02
  kord_tm = -9
  kord_mt = 9
  kord_wz = 9
  kord_tr = 9
  hydrostatic = $hydrostatic
  phys_hydrostatic = $phys_hydrostatic
  use_hydro_pressure = $use_hydro_pressure
  beta = 0.
  a_imp = 1.
  p_fac = 0.1
  k_split = $k_split
  n_split = $n_split
  nwat = 6
  na_init = $na_init
  d_ext = 0.
  dnats = 1
  fv_sg_adj = 450
  d2_bg = 0.
  nord = 2
  dddmp = 0.1
  d4_bg = 0.12
  vtdm4 = $vtdm4
  delt_max = 0.002
  ke_bg = 0.
  do_vort_damp = $do_vort_damp
  external_ic = $external_ic
  external_eta = .T.
  gfs_phil = $gfs_phil
  nggps_ic = ${nggps_ic:-".true."}
  mountain = $mountain
  ncep_ic = ${ncep_ic:-".false."}
  d_con = $d_con
  hord_mt = $hord_mt
  hord_vt = $hord_xx
  hord_tm = $hord_xx
  hord_dp = -6
  hord_tr = 8
  adjust_dry_mass = .false.
  do_sat_adj = .true.
  consv_te = $consv_te
  consv_am = .false.
  fill = .true.
  dwind_2d = .false.
  print_freq = $print_freq
  warm_start = $warm_start
  no_dycore = $no_dycore
  z_tracer = .true.
  read_increment = .false.
  res_latlon_dynamics = "fv3.increment.nc"
/
##  res_latlon_dynamics = $res_latlon_dynamics    ###DA
##  read_increment = $read_increment              ###DA

&external_ic_nml
  filtered_terrain = $filtered_terrain
  levp = $LEVS
  gfs_dwinds = $gfs_dwinds
  checker_tr = .false.
  nt_checker = 0
/

##  ntoz        = ${ntoz:-2}
##  ntcw        = ${ntcw:-3}
&gfs_physics_nml
  fhzero      = $FHZER
  ldiag3d     = ${ldiag3d:-.false.}
  fhcyc       = $FHCYC
  nst_anl     = ${nst_anl:-".true."}
  use_ufo     = ${use_ufo:-".true."}
  pre_rad     = ${pre_rad:-".false."}
  ncld        = ${ncld:-1}
  imp_physics = 11
  pdfcld      = ${pdfcld:-".false."}
  fhswr       = ${fhswr:-3600.}
  fhlwr       = ${fhlwr:-3600.}
  ialb        = ${ialb:-1}
  iems        = ${iems:-1}
  IAER        = ${iaer:-111}
  ico2        = ${ico2:-2}
  isubc_sw    = ${isubc_sw:-2}
  isubc_lw    = ${isubc_lw:-2}
  isol        = ${isol:-2}
  lwhtr       = ${lwhtr:-.true.}
  swhtr       = ${swhtr:-.true.}
  cnvgwd      = ${cnvgwd:-.true.}
  shal_cnv    = ${shal_cnv:-.true.}
  cal_pre     = ${cal_pre:-.true.}
  redrag      = ${redrag:-.true.}
  dspheat     = ${dspheat:-.true.}
  hybedmf     = ${hybedmf:-.true.}
  satmedmf    = .false.
  random_clds = ${random_clds:-.true.}
  trans_trac  = .true.
  cnvcld      = ${cnvcld:-.true.}
  imfshalcnv  = ${imfshalcnv:-2}
  imfdeepcnv  = ${imfdeepcnv:-2}
  cdmbgwd     = ${cdmbgwd:-"3.5,0.25"}
  prslrd0     = ${prslrd0:-0.}
  ivegsrc     = ${ivegsrc:-1}
  isot        = ${isot:-1}
  debug       = ${gfs_phys_debug:-".false."}
  nstf_name   = $nstf_name
  lgocart     = ${lgocart:-".true."}
  cplflx      = ${cplflx:-.false.}
  cplchm         = .true.
  ras         = ${ras:-".false"}
  iau_delthrs    = 6 
  iaufhrs        = 30
  iau_inc_files  = ' '
  psautco      = ${psautco:-"0.0008,0.0005"}
  prautco      = ${prautco:-"0.00015,0.00015"}
  lgfdlmprad   = ${lgfdlmprad:-".false."}
  effr_in      = ${effr_in:-".false."}
  fscav_aero   = "sulf:0.3","bc1:0.1","bc2:0.2","oc1:0.1","oc2:0.2",
  $gfs_physics_nml
/

&gfdl_cloud_microphysics_nml
  sedi_transport = .true.
  do_sedi_heat = .false.
  rad_snow = .true.
  rad_graupel = .true.
  rad_rain = .true.
  const_vi = .F.
  const_vs = .F.
  const_vg = .F.
  const_vr = .F.
  vi_max = 1.
  vs_max = 2.
  vg_max = 12.
  vr_max = 12.
  qi_lim = 1.
  prog_ccn = .false.
  do_qa = .true.
  fast_sat_adj = .true.
  tau_l2v = 225.
  tau_v2l = 150.
  tau_g2v = 900.
  rthresh = 10.e-6  ! This is a key parameter for cloud water
  dw_land  = 0.16
  dw_ocean = 0.10
  ql_gen = 1.0e-3
  ql_mlt = 1.0e-3
  qi0_crt = 8.0E-5
  qs0_crt = 1.0e-3
  tau_i2s = 1000.
  c_psaci = 0.05
  c_pgacs = 0.01
  rh_inc = 0.30
  rh_inr = 0.30
  rh_ins = 0.30
  ccn_l = 300.
  ccn_o = 100.
  c_paut = 0.5
  c_cracw = 0.8
  use_ppm = .false.
  use_ccn = .true.
  mono_prof = .true.
  z_slope_liq  = .true.
  z_slope_ice  = .true.
  de_ice = .false.
  fix_negative = .true.
  icloud_f = 1
  mp_time = 150.
  $gfdl_cloud_microphysics_nml
/

&interpolator_nml
  interp_method = 'conserve_great_circle'
/

&namsfc
  FNGLAC   = '${FNGLAC}' 
  FNMXIC   = '${FNMXIC}'
  FNTSFC   = '${FNTSFC}' 
  FNSNOC   = '${FNSNOC}' 
  FNZORC   = '${FNZORC}' 
  FNALBC   = '${FNALBC}' 
  FNALBC2  = '${FNALBC2}'
  FNAISC   = '${FNAISC}' 
  FNTG3C   = '${FNTG3C}' 
  FNVEGC   = '${FNVEGC}' 
  FNVETC   = '${FNVETC}' 
  FNSOTC   = '${FNSOTC}' 
  FNSMCC   = '${FNSMCC}' 
  FNMSKH   = '${FNMSKH}' 
  FNTSFA   = '${FNTSFA}' 
  FNACNA   = '${FNACNA}' 
  FNSNOA   = '${FNSNOA}' 
  FNVMNC   = '${FNVMNC}' 
  FNVMXC   = '${FNVMXC}'
  FNSLPC   = '${FNSLPC}'
  FNABSC   = '${FNABSC}'
  LDEBUG = .false.
  FSMCL(2) = 99999
  FSMCL(3) = 99999
  FSMCL(4) = 99999
  FTSFS = 90
  FAISS = 99999
  FSNOL = 99999
  FSICL = 99999
  FTSFL = 99999
  FAISL = 99999
  FVETL = 99999
  FSOTL = 99999
  FvmnL = 99999
  FvmxL = 99999
  FSLPL = 99999
  FABSL = 99999
  FSNOS = 99999
  FSICS = 99999
/
&raqmschem_nml
  chem_hist_outname = "chem_out_"
  emi_inname  = "/ships19/aqda/lenzen/RAQMSEMIS/C192/"
  fireemi_inname  = "GBBEPx/"
  emi_outname = "$emi_outname"
  gsi_path = "$gsi_path"
  wetdep_ls_opt=3
  chem_conv_tr=0
  plumerise_flag=2
  plumerisefire_frq=60
/
&chem_nml
  aer_bc_opt=1
  aer_ic_opt=1
  aer_ra_feedback=0
  aerchem_onoff=1
  bio_emiss_opt=0
  biomass_burn_opt=1
  chem_conv_tr=0
  chem_in_opt=1
  chem_opt=300
  chemdt=3
  cldchem_onoff=0
  dmsemis_opt=1
  dust_opt=5
  dust_alpha=2.3
  dust_gamma=1.0
  dust_calcdrag=1
  emiss_inpt_opt=1
  emiss_opt=5
  gas_bc_opt=1
  gas_ic_opt=1
  gaschem_onoff=1
  kemit=1
  phot_opt=1
  photdt=60
  plumerise_flag=2
  plumerisefire_frq=60
  seas_opt=2
  seas_emis_scheme=-1
  vertmix_onoff=1
  gfdlmp_onoff=1
  archive_step = -1 
  chem_hist_outname = "chem_out_"
  emi_inname  = "/ships19/aqda/lenzen/GSDCHEM_input_data/emi_C192/$SMONTH/"
  dust_inname  = "/ships19/aqda/lenzen/GSDCHEM_input_data/emi_C192/fengsha/$SMONTH/"
  fireemi_inname  = "GBBEPx/"
  emi_outname = "$emi_outname"
/
 &MOM_input_nml
         output_directory = 'OUTPUT/',
         input_filename = 'n'
         restart_input_dir = 'INPUT/',
         restart_output_dir = 'RESTART/',
         parameter_filename = 'MOM_input',
                              'MOM_override' /

 &ocean_domains_nml
 /
&nam_stochy
  lon_s=768,
  lat_s=384, 
  ntrunc=382,
  SKEBNORM=1,
  SKEB_NPASS=30,
  SKEB_VDOF=5,
  SKEB=-999,
  SKEB_TAU=2.16E4,
  SKEB_LSCALE=1000.E3,
  SHUM=-999,
  SHUM_TAU=21600,
  SHUM_LSCALE=500000,
  SPPT=-999,
  SPPT_TAU=21600,
  SPPT_LSCALE=500000,
  SPPT_LOGIT=.TRUE.,
  SPPT_SFCLIMIT=.TRUE.,
  ISEED_SHUM=1,
  ISEED_SKEB=2,
  ISEED_SPPT=3,
/

&nam_sfcperts
  NSFCPERT=6,
  PERTZ0=-999.,
  PERTSHC=-999.,
  PERTZT=-999.,
  PERTLAI=-999.,
  PERTVEGF=-999.,
  PERTALB=-999.,
  SFC_TAU=21600,
  SFC_LSCALE=500000,
  ISEED_SFC=0,
  SPPT_LAND=.FALSE.,
/
EOF



#------------------------------------------------------------------
# run the executable
cd $DATA
#export I_MPI_DEBUG=6
$NCP $FCSTEXECDIR/$FCSTEXEC $DATA/.
#$FCST_LAUNCHER ./$FCSTEXEC 1>&1 2>&2
#ajl inside
echo ajl do exec >/home/lenzen/outdoexec
#module list > /home/lenzen/module-list-before-fcst_launcher.txt # JS

#echo "$FCST_LAUNCHER ./$FCSTEXEC > $BASE_OUT/lenzen/outraqms.both.intel.submit.$PTYPE.$MODE.n$NTASKS.088.T$nth_f 2> $BASE_OUT/lenzen/errraqms.both.intel.submit.$PTYPE.$MODE.n$NTASKS.088.t$nth_f" # JS

#unset I_MPI_PMI_LIBRARY
#env >$HOME/outenvrun
#ulimit -a >$HOME/outulimit.submit2
#echo FCST_LAUNCHER $FCST_LAUNCHER >& $HOME/outlaunch
$FCST_LAUNCHER ./$FCSTEXEC > $BASE_OUT/OUTDIR/outraqms.$CASE.$CDATE.transtrac.n$NTASKS.088.T$nth_f.$FHMAX 2> $BASE_OUT/OUTDIR/errraqms.$CASE.$CDATE.transtrac.n$NTASKS.088.t$nth_f.$FHMAX

export ERR=$?
export err=$ERR
# ajl see if can keep going when freezes 3/4/2021
#$ERRSCRIPT || exit 2

#------------------------------------------------------------------
if [ $SEND = "YES" ]; then
  # Copy model output files
  cd $DATA
  for n in `seq 1 $ntiles`; do
    for file in $DATA/*.tile${n}.nc; do
#      $NCP $file $MEMDIR/.
      mv $file $MEMDIR/.
    done
  done

  # Copy model restart files
  cd $DATA/RESTART
  if [ $restart_secs -gt 0 ]; then
    restart_hrs=`echo "$restart_secs / 3600" | bc`
    RDATE=`$NDATE +$restart_hrs $CDATE`
    ryyyymmdd=`echo $RDATE | cut -c1-8`
    rhh=`echo       $RDATE | cut -c9-10`
    RMEMDIR=$ROTDIR/${PREINP}.$ryyyymmdd/$rhh/$MEMCHAR
    mkdir -p $RMEMDIR/RESTART
    for file in ${ryyyymmdd}.${rhh}0000.* ; do
      file2=$(echo $(basename $file) | sed -e "s/${ryyyymmdd}.${rhh}0000.//g")
#      $NCP $file $RMEMDIR/RESTART/$file2
      mv $file $RMEMDIR/RESTART/$file2
    done
  else
    mkdir -p $MEMDIR/RESTART
    for file in * ; do
#      $NCP $file $MEMDIR/RESTART/$file
      mv $file $MEMDIR/RESTART/$file
    done
  fi
fi

#------------------------------------------------------------------
# Clean up before leaving
if [ ${KEEPDATA:-YES} = NO ]; then rm -rf $DATA; fi

#------------------------------------------------------------------
set +x
if [ "$VERBOSE" = "YES" ] ; then
  echo $(date) EXITING $0 with return code $err >&2
fi
echo end of ajl exglobal $err
exit $err

#!/bin/ksh -x
################################################################################
# UNIX Script Documentation Block
# Script name:         exglobal_fcst_nemsfv3gfs.sh.ecf
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
# 2017-05-24  Rahul Mahajan  Updated for cycling with NEMS FV3GFS
# 2017-09-13  Fanglin Yang   Updated for using GFDL MP and Write Component
# 2019-03-21  Fanglin Yang   Add restart capability for running gfs fcst from a break point.
#
# $Id$
#
# Attributes:
#   Language: Portable Operating System Interface (POSIX) Shell
#   Machine: WCOSS-CRAY, Theia
################################################################################
#source /home/lenzen/EMC_FV3/v91/EMC_FV3GFS-GSDCHEM/release/v0/scripts/config.base
export NDATE=/data/prod/ncep_libs/intel/18.0.3/util/exec/ndate
echo CASE $CASE
echo ROTDIR $ROTDIR
echo PSLOT $PSLOT
echo top IAER $IAER
#  Set environment.
VERBOSE=${VERBOSE:-"YES"}
if [ $VERBOSE = "YES" ] ; then
  echo $(date) EXECUTING $0 $* >&2
  set -x
fi
# ajl added lgocart

machine=${machine:-"WCOSS_C"}
machine=$(echo $machine | tr '[a-z]' '[A-Z]')

# Cycling and forecast hour specific parameters
CASE=${CASE:-C96}
CDATE=${CDATE:-2016100300}
CDUMP=${CDUMP:-gdas}
FHMIN=${FHMIN:-0}
FHMAX=${FHMAX:-9}
FHOUT=${FHOUT:-6}
FHZER=${FHZER:-6}
FHCYC=${FHCYC:-24}
FHMAX_HF=${FHMAX_HF:-0}
FHOUT_HF=${FHOUT_HF:-6}
NSOUT=${NSOUT:-"-1"}
FDIAG=$FHOUT
if [ $FHMAX_HF -gt 0 -a $FHOUT_HF -gt 0 ]; then FDIAG=$FHOUT_HF; fi

PDY=$(echo $CDATE | cut -c1-8)
cyc=$(echo $CDATE | cut -c9-10)

# Directories.
pwd=$(pwd)
NWPROD=${NWPROD:-${NWROOT:-$pwd}}
HOMEgfs=${HOMEgfs:-$NWPROD}
FIX_DIR=${FIX_DIR:-$HOMEgfs/fix}
FIX_AM=${FIX_AM:-$FIX_DIR/fix_am}
FIXfv3=${FIXfv3:-$FIX_DIR/fix_fv3_gmted2010}
DATA=${DATA:-$pwd/fv3tmp$$}    # temporary running directory
mkdir -p $DATA
ROTDIR=${ROTDIR:-$pwd}         # rotating archive directory
ICSDIR=${ICSDIR:-$pwd}         # cold start initial conditions
DMPDIR=${DMPDIR:-$pwd}         # global dumps for seaice, snow and sst analysis
EMIDIR=${EMIDIR:-$pwd}         # anthro. emission 
EMITYPE=2 # we want GBBEPx
EMITYPE=${EMITYPE:-1}         # 1:MODIS, 2:GBBEPx 
# Model resolution specific parameters
DELTIM=${DELTIM:-225}
layout_x=${layout_x:-8}
layout_y=${layout_y:-16}
LEVS=${LEVS:-65}

 OUTTIME=$(( $FHOUT*3600/ $DELTIM ))
 
if [ $imp_physics -eq 99 ]; then NTRACER=0; fi
if [ $imp_physics -eq 11 ]; then NTRACER=1; fi

# Utilities
NCP=${NCP:-"/bin/cp -p"}
NLN=${NLN:-"/bin/ln -sf"}
NMV=${NMV:-"/bin/mv"}
SEND=${SEND:-"YES"}   #move final result to rotating directory
ERRSCRIPT=${ERRSCRIPT:-'eval [[ $err = 0 ]]'}
KEEPDATA=${KEEPDATA:-"NO"}

# Other options
MEMBER=${MEMBER:-"-1"} # -1: control, 0: ensemble mean, >0: ensemble member $MEMBER
ENS_NUM=${ENS_NUM:-1}  # Single executable runs multiple members (e.g. GEFS)

# Model specific stuff
FCSTEXECDIR=${FCSTEXECDIR:-$HOMEgfs/sorc/fv3gfs.fd/NEMS/exe}
FCSTEXEC=${FCSTEXEC:-fv3_gfs.x}
PARM_FV3DIAG=${PARM_FV3DIAG:-$HOMEgfs/parm/parm_fv3diag}

# Model config options
APRUN_FV3=${APRUN_FV3:-${APRUN_FCST:-${APRUN:-""}}}
NTHREADS_FV3=${NTHREADS_FV3:-${NTHREADS_FCST:-${nth_fv3:-1}}}
cores_per_node=${cores_per_node:-${npe_node_max:-24}}
ntiles=${ntiles:-6}
NTASKS_FV3=${NTASKS_FV3:-$npe_fv3}

TYPE=${TYPE:-"nh"}                  # choices:  nh, hydro
MONO=${MONO:-"non-mono"}            # choices:  mono, non-mono

CPL=${CPL:-".true."}
QUILTING=${QUILTING:-".true."}
OUTPUT_GRID=${OUTPUT_GRID:-"gaussian_grid"}
OUTPUT_FILE=${OUTPUT_FILE:-"nemsio"}
WRITE_NEMSIOFLIP=${WRITE_NEMSIOFLIP:-".true."}
WRITE_FSYNCFLAG=${WRITE_FSYNCFLAG:-".true."}

rCDUMP=${rCDUMP:-$CDUMP}

#------------------------------------------------------------------
# setup the runtime environment
if [ $machine = "WCOSS_C" ] ; then
  HUGEPAGES=${HUGEPAGES:-hugepages4M}
  . $MODULESHOME/init/sh 2>/dev/null
  module load iobuf craype-$HUGEPAGES 2>/dev/null
  export MPICH_GNI_COLL_OPT_OFF=${MPICH_GNI_COLL_OPT_OFF:-MPI_Alltoallv}
  export MKL_CBWR=AVX2
  export WRTIOBUF=${WRTIOBUF:-"4M"}
  export NC_BLKSZ=${NC_BLKSZ:-"4M"}
  export IOBUF_PARAMS="*nemsio:verbose:size=${WRTIOBUF},*:verbose:size=${NC_BLKSZ}"
fi

#-------------------------------------------------------
if [ ! -d $ROTDIR ]; then mkdir -p $ROTDIR; fi
mkdata=NO
if [ ! -d $DATA ]; then
#   mkdata=YES !lzhang
   mkdir -p $DATA
fi
# Stage the FV3 initial conditions to ROTDIR  
export OUTDIR="$ICSDIR/$CDATE/$CDUMP/$CASE/INPUT" #lzhang
COMOUT="$ROTDIR/$CDUMP.$PDY/$cyc"
[[ ! -d $COMOUT ]] && mkdir -p $COMOUT
cd $COMOUT || exit 99
[[ ! -d $INPUT ]] && $NLN $OUTDIR .
#----------------------------------------------------
cd $DATA || exit 8

mkdir -p $DATA/RESTART $DATA/INPUT
mkdir -p $DATA/GBBEPx


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
#lzhang we do not need so much restart file in the output
#if [ $CDUMP = "gfs" -a $restart_interval -gt 0 ]; then
#    RSTDIR_TMP=${RSTDIR:-$ROTDIR}/${CDUMP}.${PDY}/${cyc}/RERUN_RESTART
#    if [ ! -d $RSTDIR_TMP ]; then mkdir -p $RSTDIR_TMP ; fi
#    $NLN $RSTDIR_TMP RESTART
#else
#    mkdir -p $DATA/RESTART
#fi

#-------------------------------------------------------
# determine if restart IC exists to continue from a previous forecast
RERUN="NO"
filecount=$(find $RSTDIR_TMP -type f | wc -l) 
if [ $CDUMP = "gfs" -a $restart_interval -gt 0 -a $FHMAX -gt $restart_interval -a $filecount -gt 10 ]; then
    SDATE=$($NDATE +$FHMAX $CDATE)
    EDATE=$($NDATE +$restart_interval $CDATE)
    while [ $SDATE -gt $EDATE ]; do
        PDYS=$(echo $SDATE | cut -c1-8)
        cycs=$(echo $SDATE | cut -c9-10)
        flag1=$RSTDIR_TMP/${PDYS}.${cycs}0000.coupler.res
        flag2=$RSTDIR_TMP/coupler.res
        if [ -s $flag1 ]; then
            mv $flag1 ${flag1}.old
            if [ -s $flag2 ]; then mv $flag2 ${flag2}.old ;fi
            RERUN="YES"
            CDATE_RST=$($NDATE -$restart_interval $SDATE)
            break
        fi 
        SDATE=$($NDATE -$restart_interval $SDATE)
    done
fi

#-------------------------------------------------------
# member directory
if [ $MEMBER -lt 0 ]; then
  prefix=$CDUMP
  rprefix=$rCDUMP
  memchar=""
else
  prefix=enkf$CDUMP
  rprefix=enkf$rCDUMP
  memchar=mem$(printf %03i $MEMBER)
fi
memdir=$ROTDIR/${prefix}.$PDY/$cyc/$memchar
if [ ! -d $memdir ]; then mkdir -p $memdir; fi
assim_freq=0
GDATE=$($NDATE -$assim_freq $CDATE)
gPDY=$(echo $GDATE | cut -c1-8)
gcyc=$(echo $GDATE | cut -c9-10)
#gmemdir=$ROTDIR/${rprefix}.$gPDY/$gcyc/$memchar #lzhang
gmemdir=$ROTDIR/${prefix}.$gPDY/$gcyc/$memchar

#-------------------------------------------------------
# initial conditions
warm_start=${warm_start:-".false."}
read_increment=${read_increment:-".false."}
restart_interval=${restart_interval:-0}

# Determine if this is a warm start or cold start
if [ -f $gmemdir/RESTART/${PDY}.${cyc}0000.coupler.res ]; then
  export warm_start=".true."
fi

#-------------------------------------------------------
if [ $warm_start = ".true." -o $RERUN = "YES" ]; then
   CHEMIN=1 #lzhang
#-------------------------------------------------------
#.............................
  if [ $RERUN = "NO" ]; then
#.............................

  # Link all (except sfc_data) restart files from $gmemdir
  for file in $gmemdir/RESTART/${PDY}.${cyc}0000.*.nc; do
    file2=$(echo $(basename $file))
    file2=$(echo $file2 | cut -d. -f3-) # remove the date from file
    fsuf=$(echo $file2 | cut -d. -f1)
    if [ $fsuf != "sfc_data" ]; then
       $NLN $file $DATA/INPUT/$file2
    fi
  done

  # Link sfcanl_data restart files from $memdir
  for file in $memdir/RESTART/${PDY}.${cyc}0000.*.nc; do
    file2=$(echo $(basename $file))
    file2=$(echo $file2 | cut -d. -f3-) # remove the date from file
    fsufanl=$(echo $file2 | cut -d. -f1)
    if [ $fsufanl = "sfcanl_data" ]; then
      file2=$(echo $file2 | sed -e "s/sfcanl_data/sfc_data/g")
      $NLN $file $DATA/INPUT/$file2
    fi
  done
  
#lzhang cp atminc.nc to output dir
   cd $DATA
   $NCP ../calcinc/atminc.nc $memdir/${CDUMP}.t${cyc}z.atminc.nc

  # Handle coupler.res file for DA cycling
  if [ ${USE_COUPLER_RES:-"NO"} = "YES" ]; then
    # In DA, this is not really a "true restart",
    # and the model start time is the analysis time
    # The alternative is to replace
    # model start time with current model time in coupler.res
    file=$gmemdir/RESTART/${PDY}.${cyc}0000.coupler.res
    file2=$(echo $(basename $file))
    file2=$(echo $file2 | cut -d. -f3-) # remove the date from file
    $NLN $file $DATA/INPUT/$file2
  fi

  increment_file=$memdir/${CDUMP}.t${cyc}z.atminc.nc #lzhang 
  #increment_file=$gmemdir/${CDUMP}.t${cyc}z.atminc.nc
  if [ -f $increment_file ]; then
    $NLN $increment_file $DATA/INPUT/fv3_increment.nc
    read_increment=".true."
    res_latlon_dynamics="fv3_increment.nc"
  else
    read_increment=".false."
    res_latlon_dynamics="''"
  fi
#lzhang: link sfc data from gfs initial conditions
  for file in $memdir/INPUT/*.nc; do
    file2=$(echo $(basename $file))
    fsuf=$(echo $file2 | cut -c1-3)
    if [ $fsuf = $fsuf = "sfc" ]; then
      $NLN $file $DATA/INPUT/$file2
    fi
  done
#lzhang
#.............................
  else  ##RERUN                         

    PDYT=$(echo $CDATE_RST | cut -c1-8)
    cyct=$(echo $CDATE_RST | cut -c9-10)
    for file in $RSTDIR_TMP/${PDYT}.${cyct}0000.*; do
      file2=$(echo $(basename $file))
      file2=$(echo $file2 | cut -d. -f3-)
      $NLN $file $DATA/INPUT/$file2
    done

  fi
#.............................
else ## cold start                            
  CHEMIN=0 #lzhang
  echo "COLD START"
# ajl need to change to our path on s4

 if [ -d $IC_DIR/${CASE}_$CDATE ]; then
  $NCP $IC_DIR/${CASE}_$CDATE/* $DATA/INPUT/.
 else # ajl
   for file in $memdir/INPUT/*.nc; do
     file2=$(echo $(basename $file))
     fsuf=$(echo $file2 | cut -c1-3)
     if [ $fsuf = "gfs" -o $fsuf = "sfc" ]; then
       $NLN $file $DATA/INPUT/$file2
     fi
    done
  fi # ajl
  CHEMIN=1 # ajl for now since we put it all in with met

#-------------------------------------------------------
fi 
#-------------------------------------------------------


nfiles=$(ls -1 $DATA/INPUT/* | wc -l)
if [ $nfiles -le 0 ]; then
  echo "Initial conditions must exist in $DATA/INPUT, ABORT!"
  msg=”"Initial conditions must exist in $DATA/INPUT, ABORT!"
  postmsg "$jlogfile" "$msg"
  exit 1
fi

#--------------------------------------------------------------------------
# Grid and orography data
for n in $(seq 1 $ntiles); do
  $NLN $FIXfv3/$CASE/${CASE}_grid.tile${n}.nc     $DATA/INPUT/${CASE}_grid.tile${n}.nc
  $NLN $FIXfv3/$CASE/${CASE}_oro_data.tile${n}.nc $DATA/INPUT/oro_data.tile${n}.nc
done
$NLN $FIXfv3/$CASE/${CASE}_mosaic.nc  $DATA/INPUT/grid_spec.nc

# GFS standard input data
echo IAER $IAER
export IAER=111
IALB=${IALB:-1}
IEMS=${IEMS:-1}
ISOL=${ISOL:-2}
IAER=${IAER:-111}
ICO2=${ICO2:-2}

if [ ${new_o3forc:-YES} = YES ]; then
    O3FORC=ozprdlos_2015_new_sbuvO3_tclm15_nuchem.f77
else
    O3FORC=global_o3prdlos.f77
fi
H2OFORC=${H2OFORC:-"global_h2o_pltc.f77"}
$NLN $FIX_AM/${O3FORC}                         $DATA/global_o3prdlos.f77
$NLN $FIX_AM/${H2OFORC}                        $DATA/global_h2oprdlos.f77
$NLN $FIX_AM/global_solarconstant_noaa_an.txt  $DATA/solarconstant_noaa_an.txt
$NLN $FIX_AM/global_sfc_emissivity_idx.txt     $DATA/sfc_emissivity_idx.txt

$NLN $FIX_AM/global_co2historicaldata_glob.txt $DATA/co2historicaldata_glob.txt
$NLN $FIX_AM/co2monthlycyc.txt                 $DATA/co2monthlycyc.txt
if [ $ICO2 -gt 0 ]; then
  for file in $(ls $FIX_AM/fix_co2_proj/global_co2historicaldata*) ; do
    $NLN $file $DATA/$(echo $(basename $file) | sed -e "s/global_//g")
  done
fi

$NLN $FIX_AM/global_climaeropac_global.txt     $DATA/aerosol.dat
if [ $IAER -gt 0 ] ; then
  for file in $(ls $FIX_AM/global_volcanic_aerosols*) ; do
    $NLN $file $DATA/$(echo $(basename $file) | sed -e "s/global_//g")
  done
fi
#------------------------------------------------------------------
# changeable parameters
# dycore definitions
res=$(echo $CASE |cut -c2-5)
resp=$((res+1))
npx=$resp
npy=$resp
npz=$((LEVS-1))
io_layout="1,1"
#ncols=$(( (${npx}-1)*(${npy}-1)*3/2 ))

# spectral truncation and regular grid resolution based on FV3 resolution
JCAP_CASE=$((2*res-2))
LONB_CASE=$((4*res))
LATB_CASE=$((2*res))

JCAP=${JCAP:-$JCAP_CASE}
LONB=${LONB:-$LONB_CASE}
LATB=${LATB:-$LATB_CASE}

LONB_IMO=${LONB_IMO:-$LONB_CASE}
LATB_JMO=${LATB_JMO:-$LATB_CASE}

# Fix files
FNGLAC=${FNGLAC:-"$FIX_AM/global_glacier.2x2.grb"}
FNMXIC=${FNMXIC:-"$FIX_AM/global_maxice.2x2.grb"}
FNTSFC=${FNTSFC:-"$FIX_AM/RTGSST.1982.2012.monthly.clim.grb"}
FNSNOC=${FNSNOC:-"$FIX_AM/global_snoclim.1.875.grb"}
FNZORC=${FNZORC:-"igbp"}
FNALBC2=${FNALBC2:-"$FIX_AM/global_albedo4.1x1.grb"}
FNAISC=${FNAISC:-"$FIX_AM/CFSR.SEAICE.1982.2012.monthly.clim.grb"}
FNTG3C=${FNTG3C:-"$FIX_AM/global_tg3clim.2.6x1.5.grb"}
FNVEGC=${FNVEGC:-"$FIX_AM/global_vegfrac.0.144.decpercent.grb"}
FNMSKH=${FNMSKH:-"$FIX_AM/seaice_newland.grb"}
FNVMNC=${FNVMNC:-"$FIX_AM/global_shdmin.0.144x0.144.grb"}
FNVMXC=${FNVMXC:-"$FIX_AM/global_shdmax.0.144x0.144.grb"}
FNSLPC=${FNSLPC:-"$FIX_AM/global_slope.1x1.grb"}
FNALBC=${FNALBC:-"$FIX_AM/global_snowfree_albedo.bosu.t${JCAP}.${LONB}.${LATB}.rg.grb"}
FNVETC=${FNVETC:-"$FIX_AM/global_vegtype.igbp.t${JCAP}.${LONB}.${LATB}.rg.grb"}
FNSOTC=${FNSOTC:-"$FIX_AM/global_soiltype.statsgo.t${JCAP}.${LONB}.${LATB}.rg.grb"}
FNABSC=${FNABSC:-"$FIX_AM/global_mxsnoalb.uariz.t${JCAP}.${LONB}.${LATB}.rg.grb"}
FNSMCC=${FNSMCC:-"$FIX_AM/global_soilmgldas.statsgo.t${JCAP}.${LONB}.${LATB}.grb"}

# If the appropriate resolution fix file is not present, use the highest resolution available (T1534)
[[ ! -f $FNALBC ]] && FNALBC="$FIX_AM/global_snowfree_albedo.bosu.t1534.3072.1536.rg.grb"
[[ ! -f $FNVETC ]] && FNVETC="$FIX_AM/global_vegtype.igbp.t1534.3072.1536.rg.grb"
[[ ! -f $FNSOTC ]] && FNSOTC="$FIX_AM/global_soiltype.statsgo.t1534.3072.1536.rg.grb"
[[ ! -f $FNABSC ]] && FNABSC="$FIX_AM/global_mxsnoalb.uariz.t1534.3072.1536.rg.grb"
[[ ! -f $FNSMCC ]] && FNSMCC="$FIX_AM/global_soilmgldas.statsgo.t1534.3072.1536.grb"

# NSST Options
# nstf_name contains the NSST related parameters
# nstf_name(1) : NST_MODEL (NSST Model) : 0 = OFF, 1 = ON but uncoupled, 2 = ON and coupled
# nstf_name(2) : NST_SPINUP : 0 = OFF, 1 = ON,
# nstf_name(3) : NST_RESV (Reserved, NSST Analysis) : 0 = OFF, 1 = ON
# nstf_name(4) : ZSEA1 (in mm) : 0
# nstf_name(5) : ZSEA2 (in mm) : 0
# nst_anl      : .true. or .false., NSST analysis over lake
NST_MODEL=${NST_MODEL:-0}
NST_SPINUP=${NST_SPINUP:-0}
NST_RESV=${NST_RESV-0}
ZSEA1=${ZSEA1:-0}
ZSEA2=${ZSEA2:-0}
nstf_name=${nstf_name:-"$NST_MODEL,$NST_SPINUP,$NST_RESV,$ZSEA1,$ZSEA2"}
nst_anl=${nst_anl:-".false."}


# blocking factor used for threading and general physics performance
#nyblocks=`expr \( $npy - 1 \) \/ $layout_y `
#nxblocks=`expr \( $npx - 1 \) \/ $layout_x \/ 32`
#if [ $nxblocks -le 0 ]; then nxblocks=1 ; fi
blocksize=${blocksize:-32}

# the pre-conditioning of the solution
# =0 implies no pre-conditioning
# >0 means new adiabatic pre-conditioning
# <0 means older adiabatic pre-conditioning
na_init=${na_init:-1}
[[ $warm_start = ".true." ]] && na_init=0

# variables for controlling initialization of NCEP/NGGPS ICs
filtered_terrain=${filtered_terrain:-".true."}
gfs_dwinds=${gfs_dwinds:-".true."}

# various debug options
no_dycore=${no_dycore:-".false."}
dycore_only=${adiabatic:-".false."}
chksum_debug=${chksum_debug:-".false."}
print_freq=${print_freq:-6}

if [ ${TYPE} = "nh" ]; then # non-hydrostatic options

  hydrostatic=".false."
  phys_hydrostatic=".false."     # enable heating in hydrostatic balance in non-hydrostatic simulation
  use_hydro_pressure=".false."   # use hydrostatic pressure for physics
  if [ $warm_start = ".true." ]; then
    make_nh=".false."              # restarts contain non-hydrostatic state
  else
    make_nh=".true."               # re-initialize non-hydrostatic state
  fi

else # hydrostatic options

  hydrostatic=".true."
  phys_hydrostatic=".false."     # ignored when hydrostatic = T
  use_hydro_pressure=".false."   # ignored when hydrostatic = T
  make_nh=".false."              # running in hydrostatic mode

fi

# Conserve total energy as heat globally
consv_te=${consv_te:-1.} # range 0.-1., 1. will restore energy to orig. val. before physics

# time step parameters in FV3
k_split=${k_split:-2}
n_split=${n_split:-6}

if [ $(echo $MONO | cut -c-4) = "mono" ];  then # monotonic options

  d_con=${d_con_mono:-"0."}
  do_vort_damp=".false."
  if [ ${TYPE} = "nh" ]; then # non-hydrostatic
    hord_mt=${hord_mt_nh_mono:-"10"}
    hord_xx=${hord_xx_nh_mono:-"10"}
  else # hydrostatic
    hord_mt=${hord_mt_hydro_mono:-"10"}
    hord_xx=${hord_xx_hydro_mono:-"10"}
  fi

else # non-monotonic options

  d_con=${d_con_nonmono:-"1."}
  do_vort_damp=".true."
  if [ ${TYPE} = "nh" ]; then # non-hydrostatic
    hord_mt=${hord_mt_nh_nonmono:-"5"}
    hord_xx=${hord_xx_nh_nonmono:-"5"}
  else # hydrostatic
    hord_mt=${hord_mt_hydro_nonmono:-"10"}
    hord_xx=${hord_xx_hydro_nonmono:-"10"}
  fi

fi

if [ $(echo $MONO | cut -c-4) != "mono" -a $TYPE = "nh" ]; then
  vtdm4=${vtdm4_nh_nonmono:-"0.06"}
else
  vtdm4=${vtdm4:-"0.05"}
fi

if [ $warm_start = ".true." ]; then # warm start from restart file

  nggps_ic=".false."
  ncep_ic=".false."
  external_ic=".false."
  mountain=".true."
  if [ $read_increment = ".true." ]; then # add increment on the fly to the restarts
    res_latlon_dynamics="fv3_increment.nc"
  else
    res_latlon_dynamics='""'
  fi

else # CHGRES'd GFS analyses

  nggps_ic=${nggps_ic:-".true."}
  ncep_ic=${ncep_ic:-".false."}
  external_ic=".true."
  mountain=".false."
  read_increment=".false."
  res_latlon_dynamics='""'

fi

# Stochastic Physics Options
if [ ${SET_STP_SEED:-"YES"} = "YES" ]; then
  ISEED_SKEB=$((CDATE*1000 + MEMBER*10 + 1))
  ISEED_SHUM=$((CDATE*1000 + MEMBER*10 + 2))
  ISEED_SPPT=$((CDATE*1000 + MEMBER*10 + 3))
else
  ISEED=${ISEED:-0}
fi
DO_SKEB=${DO_SKEB:-"NO"}
DO_SPPT=${DO_SPPT:-"NO"}
DO_SHUM=${DO_SHUM:-"NO"}
JCAP_STP=${JCAP_STP:-$JCAP_CASE}
LONB_STP=${LONB_STP:-$LONB_CASE}
LATB_STP=${LATB_STP:-$LATB_CASE}

# build the date for curr_date and diag_table from CDATE
SYEAR=$(echo  $CDATE | cut -c1-4)
SMONTH=$(echo $CDATE | cut -c5-6)
SDAY=$(echo   $CDATE | cut -c7-8)
SHOUR=$(echo  $CDATE | cut -c9-10)
curr_date="${SYEAR},${SMONTH},${SDAY},${SHOUR},0,0"
rsecs=$((restart_interval*3600))
restart_secs=${rsecs:-0}

# copy over the tables
DIAG_TABLE=${DIAG_TABLE:-$PARM_FV3DIAG/diag_table}
DATA_TABLE=${DATA_TABLE:-$PARM_FV3DIAG/data_table}
FIELD_TABLE=${FIELD_TABLE:-$PARM_FV3DIAG/field_table}

# build the diag_table with the experiment name and date stamp
cat > diag_table << EOF
FV3 Forecast
$SYEAR $SMONTH $SDAY $SHOUR 0 0
EOF
cat $DIAG_TABLE >> diag_table

$NCP $DATA_TABLE  data_table
$NCP $FIELD_TABLE field_table
# for raqmschem

$NCP $AODLUT/* .

#------------------------------------------------------------------
rm -f nems.configure
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
  @$DELTIM
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

rm -f model_configure
cat > model_configure <<EOF
total_member:            $ENS_NUM
print_esmf:              ${print_esmf:-.true.}
PE_MEMBER01:             $NTASKS_FV3
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
cpl:                     $CPL
memuse_verbose:          ${memuse_verbose:-".false."}
atmos_nthreads:          $NTHREADS_FV3
use_hyper_thread:        ${hyperthread:-".false."}
ncores_per_node:         $cores_per_node
restart_interval:        $restart_interval
output_1st_tstep_rst:    .false.

quilting:                $QUILTING
write_groups:            ${WRITE_GROUP:-1}
write_tasks_per_group:   ${WRTTASK_PER_GROUP:-6}
num_files:               ${NUM_FILES:-2}
filename_base:           'atm' 'sfc'
output_grid:             $OUTPUT_GRID
output_file:             $OUTPUT_FILE
write_nemsioflip:        $WRITE_NEMSIOFLIP
write_fsyncflag:         $WRITE_FSYNCFLAG
imo:                     $LONB_IMO
jmo:                     $LATB_JMO

nfhout:                  $FHOUT
nfhmax_hf:               $FHMAX_HF
nfhout_hf:               $FHOUT_HF
nsout:                   $NSOUT
EOF

#&coupler_nml
#  months = ${months:-0}
#  days = ${days:-$((FHMAX/24))}
#  hours = ${hours:-$((FHMAX-24*(FHMAX/24)))}
#  dt_atmos = $DELTIM
#  dt_ocean = $DELTIM
#  current_date = $curr_date
#  calendar = 'julian'
#  memuse_verbose = .false.
#  atmos_nthreads = $NTHREADS_FV3
#  use_hyper_thread = ${hyperthread:-".false."}
#  ncores_per_node = $cores_per_node
#  restart_secs = $restart_secs
#  $coupler_nml
#/

cat > input.nml <<EOF
&amip_interp_nml
  interp_oi_sst = .true.
  use_ncep_sst = .true.
  use_ncep_ice = .false.
  no_anom_sst = .false.
  data_set = 'reynolds_oi'
  date_out_of_range = 'climo'
  $amip_interp_nml
/

&atmos_model_nml
  blocksize = $blocksize
  chksum_debug = $chksum_debug
  dycore_only = $dycore_only
  fdiag = $FDIAG
  fhmax = $FHMAX
  fhout = $FHOUT
  fhmaxhf = $FHMAX_HF
  fhouthf = $FHOUT_HF
  $atmos_model_nml
/

&diag_manager_nml
  prepend_date = .false.
  $diag_manager_nml
/

&fms_io_nml
  checksum_required = .false.
  max_files_r = 100
  max_files_w = 100
  $fms_io_nml
/

&fms_nml
  clock_grain = 'ROUTINE'
  domains_stack_size = ${domains_stack_size:-300000000}
  print_memory_usage = ${print_memory_usage:-".false."}
  $fms_nml
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
  n_sponge = ${n_sponge:-"10"}
  nudge_qv = ${nudge_qv:-".true."}
  nudge_dz = ${nudge_dz:-".false."}
  tau = ${tau:-10.}
  rf_cutoff = ${rf_cutoff:-"7.5e2"}
  d2_bg_k1 = ${d2_bg_k1:-"0.15"}
  d2_bg_k2 = ${d2_bg_k2:-"0.02"}
  kord_tm = ${kord_tm:-"-9"}
  kord_mt = ${kord_mt:-"9"}
  kord_wz = ${kord_wz:-"9"}
  kord_tr = ${kord_tr:-"9"}
  hydrostatic = $hydrostatic
  phys_hydrostatic = $phys_hydrostatic
  use_hydro_pressure = $use_hydro_pressure
  beta = 0.
  a_imp = 1.
  p_fac = 0.1
  k_split = $k_split
  n_split = $n_split
  nwat = ${nwat:-2}
  na_init = $na_init
  d_ext = 0.
  dnats = ${dnats:-0}
  fv_sg_adj = ${fv_sg_adj:-"450"}
  d2_bg = 0.
  nord = ${nord:-3}
  dddmp = ${dddmp:-0.2}
  d4_bg = ${d4_bg:-0.15}
  vtdm4 = $vtdm4
  delt_max = ${delt_max:-"0.002"}
  ke_bg = 0.
  do_vort_damp = $do_vort_damp
  external_ic = $external_ic
  external_eta = ${external_eta:-.true.}
  gfs_phil = ${gfs_phil:-".false."}
  nggps_ic = $nggps_ic
  mountain = $mountain
  ncep_ic = $ncep_ic
  d_con = $d_con
  hord_mt = $hord_mt
  hord_vt = $hord_xx
  hord_tm = $hord_xx
  hord_dp = -$hord_xx
  hord_tr = ${hord_tr:-"8"}
  adjust_dry_mass = ${adjust_dry_mass:-".false."}
  consv_te = $consv_te
  do_sat_adj = ${do_sat_adj:-".false."}
  consv_am = .false.
  fill = .true.
  dwind_2d = .false.
  print_freq = $print_freq
  warm_start = $warm_start
  no_dycore = $no_dycore
  z_tracer = .true.
  agrid_vel_rst = ${agrid_vel_rst:-".true."}
  read_increment = $read_increment
  res_latlon_dynamics = $res_latlon_dynamics
  $fv_core_nml
/

&external_ic_nml
  filtered_terrain = $filtered_terrain
  levp = $LEVS
  gfs_dwinds = $gfs_dwinds
  checker_tr = .false.
  nt_checker = 0
  $external_ic_nml
/

&gfs_physics_nml
  fhzero       = $FHZER
  h2o_phys     = ${h2o_phys:-".true."}
  ldiag3d      = ${ldiag3d:-".false."}
  fhcyc        = $FHCYC
  use_ufo      = ${use_ufo:-".true."}
  pre_rad      = ${pre_rad:-".false."}
  ncld         = ${ncld:-1}
  imp_physics  = ${imp_physics:-"99"}
  pdfcld       = ${pdfcld:-".false."}
  fhswr        = ${FHSWR:-"3600."}
  fhlwr        = ${FHLWR:-"3600."}
  ialb         = $IALB
  iems         = $IEMS
  iaer         = $IAER
  ico2         = $ICO2
  lgocart      = ${lgocart:-".true."}
  isubc_sw     = ${isubc_sw:-"2"}
  isubc_lw     = ${isubc_lw:-"2"}
  isol         = $ISOL
  lwhtr        = ${lwhtr:-".true."}
  swhtr        = ${swhtr:-".true."}
  cnvgwd       = ${cnvgwd:-".true."}
  shal_cnv     = ${shal_cnv:-".true."}
  cal_pre      = ${cal_pre:-".true."}
  redrag       = ${redrag:-".true."}
  dspheat      = ${dspheat:-".true."}
  hybedmf      = ${hybedmf:-".true."}
  random_clds  = ${random_clds:-".true."}
  trans_trac   = ${trans_trac:-".true."}
  cnvcld       = ${cnvcld:-".true."}
  imfshalcnv   = ${imfshalcnv:-"2"}
  imfdeepcnv   = ${imfdeepcnv:-"2"}
  cdmbgwd      = ${cdmbgwd:-"3.5,0.25"}
  prslrd0      = ${prslrd0:-"0."}
  ivegsrc      = ${ivegsrc:-"1"}
  isot         = ${isot:-"1"}
  debug        = ${gfs_phys_debug:-".false."}
  nstf_name    = $nstf_name
  cplflx       = .false.
  cplchm       = .true.
  nst_anl      = $nst_anl
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
  $interpolator_nml
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
  LDEBUG = ${LDEBUG:-".false."}
  FSMCL(2) = ${FSMCL2:-99999}
  FSMCL(3) = ${FSMCL3:-99999}
  FSMCL(4) = ${FSMCL4:-99999}
  FTSFS = ${FTSFS:-90}
  FAISL = ${FAISL:-99999}
  FAISS = ${FAISS:-99999}
  FSNOL = ${FSNOL:-99999}
  FSNOS = ${FSNOS:-99999}
  FSICL = 99999
  FSICS = 99999
  FTSFL = 99999
  FVETL = 99999
  FSOTL = 99999
  FvmnL = 99999
  FvmxL = 99999
  FSLPL = 99999
  FABSL = 99999
  $namsfc_nml
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
  plumerisefire_frq=30
/

&chem_nml
  aer_bc_opt=1
  aer_ic_opt=1
  aer_ra_feedback=2
  aerchem_onoff=1
  bio_emiss_opt=0
  biomass_burn_opt=1
  chem_conv_tr=0
  chem_in_opt=$CHEMIN
  chem_opt=300
  chemdt=3
  cldchem_onoff=0
  dmsemis_opt=1
  dust_opt=5
  dust_alpha = 2.3
  dust_gamma = 1.0
  dust_calcdrag=1
  emiss_inpt_opt=1
  emiss_opt=5
  gas_bc_opt=1
  gas_ic_opt=1
  gaschem_onoff=1
  kemit=1
  phot_opt=1
  photdt=60
  plumerisefire_frq=30
  PLUMERISE_flag=$EMITYPE
  seas_opt=2
  seas_emis_scheme=-1
  seas_emis_scale=1.,1.,1.,1.,1.
  vertmix_onoff=1
  gfdlmp_onoff=$NTRACER
  archive_step =-1
  chem_hist_outname = "chem_out_"
  emi_inname  = "/ships19/aqda/lenzen/CEDS_2019/GEFS-Aerosol_emissions/emi_${CASE}_RDRAG/$SMONTH"
  dust_inname = "/ships19/aqda/lenzen/CEDS_2019/GEFS-Aerosol_emissions/emi_${CASE}_RDRAG/fengsha/$SMONTH"
  fireemi_inname  = "GBBEPx"
  emi_outname = "./"
  $chem_nml
/

&fv_grid_nml
  grid_file = 'INPUT/grid_spec.nc'
  $fv_grid_nml
/
EOF

# Add namelist for stochastic physics options
echo "" >> input.nml
if [ $MEMBER -gt 0 ]; then

    cat >> input.nml << EOF
&nam_stochy
  ntrunc = $JCAP_STP
  lon_s = $LONB_STP
  lat_s = $LATB_STP
EOF

  if [ $DO_SKEB = "YES" ]; then
    cat >> input.nml << EOF
  skeb = $SKEB
  iseed_skeb = ${ISEED_SKEB:-$ISEED}
  skeb_tau = ${SKEB_TAU:-"-999."}
  skeb_lscale = ${SKEB_LSCALE:-"-999."}
  skebnorm = ${SKEBNORM:-"1"}
EOF
  fi

  if [ $DO_SHUM = "YES" ]; then
    cat >> input.nml << EOF
  shum = $SHUM
  iseed_shum = ${ISEED_SHUM:-$ISEED}
  shum_tau = ${SHUM_TAU:-"-999."}
  shum_lscale = ${SHUM_LSCALE:-"-999."}
EOF
  fi

  if [ $DO_SPPT = "YES" ]; then
    cat >> input.nml << EOF
  sppt = $SPPT
  iseed_sppt = ${ISEED_SPPT:-$ISEED}
  sppt_tau = ${SPPT_TAU:-"-999."}
  sppt_lscale = ${SPPT_LSCALE:-"-999."}
  sppt_logit = ${SPPT_LOGIT:-".true."}
  sppt_sfclimit = ${SPPT_SFCLIMIT:-".true."}
  use_zmtnblck = ${use_zmtnblck:-".true."}
EOF
  fi

  cat >> input.nml << EOF
  $nam_stochy_nml
/
EOF


    cat >> input.nml << EOF
&nam_sfcperts
  $nam_sfcperts_nml
/
EOF

else

  cat >> input.nml << EOF
&nam_stochy
/
&nam_sfcperts
/
EOF

fi


#------------------------------------------------------------------
# make symbolic links to write forecast files directly in memdir
cd $DATA
if [ $QUILTING = ".true." -a $OUTPUT_GRID = "gaussian_grid" ]; then
  fhr=$FHMIN
  while [ $fhr -le $FHMAX ]; do
    FH3=$(printf %03i $fhr)
    atmi=atmf${FH3}.$OUTPUT_FILE
    sfci=sfcf${FH3}.$OUTPUT_FILE
    logi=logf${FH3}
    atmo=$memdir/${CDUMP}.t${cyc}z.atmf${FH3}.$OUTPUT_FILE
    sfco=$memdir/${CDUMP}.t${cyc}z.sfcf${FH3}.$OUTPUT_FILE
    logo=$memdir/${CDUMP}.t${cyc}z.logf${FH3}.$OUTPUT_FILE
    eval $NLN $atmo $atmi
    eval $NLN $sfco $sfci
    eval $NLN $logo $logi
    FHINC=$FHOUT
    if [ $FHMAX_HF -gt 0 -a $FHOUT_HF -gt 0 -a $fhr -lt $FHMAX_HF ]; then
      FHINC=$FHOUT_HF
    fi
    fhr=$((fhr+FHINC))
  done
else
  for n in $(seq 1 $ntiles); do
    eval $NLN nggps2d.tile${n}.nc       $memdir/nggps2d.tile${n}.nc
    eval $NLN nggps3d.tile${n}.nc       $memdir/nggps3d.tile${n}.nc
    eval $NLN grid_spec.tile${n}.nc     $memdir/grid_spec.tile${n}.nc
    eval $NLN atmos_static.tile${n}.nc  $memdir/atmos_static.tile${n}.nc
    eval $NLN atmos_4xdaily.tile${n}.nc $memdir/atmos_4xdaily.tile${n}.nc
  done
fi

# Copy namelist file
$NCP input.nml $memdir

#------------------------------------------------------------------
# run the executable

$NCP $FCSTEXECDIR/$FCSTEXEC $DATA/.
export OMP_NUM_THREADS=$NTHREADS_FV3
export MKL_NUM_THREADS=0
#$APRUN_FV3 /bin/env KMP_AFFINITY=scatter KMP_NUM_THREADS=1 $DATA/$FCSTEXEC 1>&1 2>&2
$APRUN_FV3  $DATA/$FCSTEXEC >$BASE_OUT/OUTDIR/outraqms.$CASE.$CDATE.transtrac.n$NTASKS.088.T$nth_f.$FHMAX 2> $BASE_OUT/OUTDIR/errraqms.$CASE.$CDATE.transtrac.n$NTASKS.088.t$nth_f.$FHMAX
export ERR=$?
export err=$ERR
$ERRSCRIPT || exit $err

#------------------------------------------------------------------
if [ $SEND = "YES" ]; then
  # Copy model restart files
  cd $DATA/RESTART
  mkdir -p $memdir/RESTART

  # Only save restarts at single time in RESTART directory
  # Either at restart_interval or at end of the forecast
  if [ $restart_interval -eq 0 -o $restart_interval -eq $FHMAX ]; then

    # Add time-stamp to restart files at FHMAX
    RDATE=$($NDATE +$FHMAX $CDATE)
    rPDY=$(echo $RDATE | cut -c1-8)
    rcyc=$(echo $RDATE | cut -c9-10)
    for file in $(ls * | grep -v 0000); do
      $NMV $file ${rPDY}.${rcyc}0000.$file
      $NCP ${rPDY}.${rcyc}0000.$file $memdir/RESTART/${rPDY}.${rcyc}0000.$file #lzhang
    done

    RSTDIR_TMP=${RSTDIR:-$ROTDIR}/${CDUMP}.${PDY}/${cyc}/RERUN_RESTART #lzhang
    #if [ ! -d $RSTDIR_TMP ]; then mkdir -p $RSTDIR_TMP ; fi #lzhang
    $NLN $memdir/RESTART  $RSTDIR_TMP

  else

    # time-stamp exists at restart_interval time, just copy
    #lzhang, copy all the restart file
    for ns in $(seq 1 $FHMAX); do  
    RDATE=$($NDATE +$ns $CDATE)
    #RDATE=$($NDATE +$restart_interval $CDATE)
    rPDY=$(echo $RDATE | cut -c1-8)
    rcyc=$(echo $RDATE | cut -c9-10)
    for file in ${rPDY}.${rcyc}0000.* ; do
      $NCP $file $memdir/RESTART/$file
    done
    done

  fi

fi

#------------------------------------------------------------------
# Clean up before leaving
if [ $mkdata = "YES" ]; then rm -rf $DATA; fi

#------------------------------------------------------------------
set +x
if [ $VERBOSE = "YES" ] ; then
  echo $(date) EXITING $0 with return code $err >&2
fi
exit 0

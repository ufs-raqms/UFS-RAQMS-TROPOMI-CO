#!/bin/bash
set -eux

source rt_utils.sh
source atparse.bash

mkdir -p ${RUNDIR}
cd $RUNDIR

###############################################################################
# Make configure and run files
###############################################################################

# FV3 executable:
#cp ${PATHRT}/$FV3X                                 fv3.exe

# modulefile for FV3 prerequisites:
#cp ${PATHRT}/modules.fv3_${COMPILE_NR}             modules.fv3

# Get the shell file that loads the "module" command and purges modules:
cp ${PATHRT}/../NEMS/src/conf/module-setup.sh.inc  module-setup.sh
echo PATHRT $PATHRT
RUND="${RUNDIR}"
export PATHTR=/home/lenzen/EMC_FV3/v91/EMC_FV3GFS-GSDCHEM/
export PATHRT=/home/lenzen/EMC_FV3/v91/EMC_FV3GFS-GSDCHEM/tests/
#SRCD="${PATHTR}"
SRCD=${PATHTR}
echo RUN

#atparse < ${PATHRT}/fv3_conf/${FV3_RUN:-fv3gsdchem_run.IN} > fv3gsdchem_run
echo INPUT

atparse < ${PATHTR}/parm/${INPUT_NML:-input.nml.IN} > input.nml
echo MODEL_CONGIFURE

atparse < ${PATHTR}/parm/${MODEL_CONFIGURE:-model_configure.IN} > model_configure
exit 0

if [[ "Q${INPUT_NEST02_NML:-}" != Q ]] ; then
    atparse < ${PATHTR}/parm/${INPUT_NEST02_NML} > input_nest02.nml
fi

source ./fv3gsdchem_run

echo "DIAG_TABLE = ${DIAG_TABLE:-}"
echo "FIELD_TABLE = ${FIELD_TABLE:-}"
echo "NEMS_CONFIGURE = ${NEMS_CONFIGURE:-}"

if [[ "Q${DIAG_TABLE:-}" != Q ]] ; then
    cp ${PATHTR}/parm/${DIAG_TABLE}              diag_table
fi

if [[ "Q${FIELD_TABLE:-}" != Q ]] ; then
    cp ${PATHTR}/parm/${FIELD_TABLE}             field_table
fi

if [[ "Q${NEMS_CONFIGURE:-}" != Q ]] ; then
    atparse < ${PATHTR}/parm/${NEMS_CONFIGURE} > nems.configure
fi

if [[ $SCHEDULER = 'moab' ]]; then
  atparse < $PATHRT/fv3_conf/fv3_msub.IN > job_card
elif [[ $SCHEDULER = 'pbs' ]]; then
  NODES=$(( TASKS / TPN ))
  if (( NODES * TPN < TASKS )); then
    NODES=$(( NODES + 1 ))
  fi
  atparse < $PATHRT/fv3_conf/fv3_qsub.IN > job_card
elif [[ $SCHEDULER = 'slurm' ]]; then
  atparse < $PATHRT/fv3_conf/fv3_slurm.IN > job_card
elif [[ $SCHEDULER = 'lsf' ]]; then
  if (( TASKS < TPN )); then
    TPN=${TASKS}
  fi
  atparse < $PATHRT/fv3_conf/fv3_bsub.IN > job_card
fi

################################################################################
# Submit test
################################################################################

if [[ $ROCOTO = 'false' ]]; then
  submit_and_wait job_card
else
  chmod u+x job_card
  ./job_card
fi

check_results

################################################################################
# End test
################################################################################

exit 0

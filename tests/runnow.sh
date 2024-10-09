#!/bin/bash
export DISKNM=mine
export RT_SUFFIX=_prod
export BL_SUFFIX=_ccpp
export LOG_DIR=mine
./run_test.sh $PWD $PWD/rundir cpld_fv3_gfdlmp_gsdchem_gbbepx_frp_fengsha newrun intel

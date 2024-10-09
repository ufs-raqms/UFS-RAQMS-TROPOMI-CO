#!/bin/bash -x
export CDATE=2919080808
export HH=`echo $CDATE | cut -c9-10`
export HH=3$HH
echo HH $HH
export HH3=$((HH % 3))
echo HH3 $HHE3

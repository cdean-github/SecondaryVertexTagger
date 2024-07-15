#!/bin/bash

source /opt/sphenix/core/bin/sphenix_setup.sh -n new

#export HOME=/sphenix/u/${LOGNAME}
export HOME=/sphenix/u/cdean
export SPHENIX=$HOME/sPHENIX
export MYINSTALL=$SPHENIX/install
export LD_LIBRARY_PATH=$MYINSTALL/lib:$LD_LIBRARY_PATH
export ROOT_INCLUDE_PATH=$MYINSTALL/include:$ROOT_INCLUDE_PATH

source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL

inputFile="/direct/sphenix+tg+tg01/hf/cdean/selected_hf_simulations/siliconOnly/D2Kpi_20240711/D2Kpi_DST_00000.root"
if [[ "$1" != "" ]]; then
    inputFile="$1"
fi

echo running: runSV.sh 
root.exe -q -b Fun4All_myTester.C\(\"$inputFile\"\)
echo Script done

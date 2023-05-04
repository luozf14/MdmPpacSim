#!/bin/sh

stage1_config="stage1.json"
stage2_config="stage2.json"

./MdmPpacSim $stage1_config $1
if [ $? -eq 0 ]; then
    echo "---> Success: Stage1"
    ./MdmPpacSim $stage2_config $1
else
    echo "---> Error: Stage1"
    exit
fi
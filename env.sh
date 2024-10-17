#!/bin/bash
# tested on linux only, might work on macos too

if [ -z "$PATHTOOLS" ]; then
    # PATHTOOLS is not defined, so check if 'vmd' exists in the system
    read -p "Enter path of tools directory: " PATHTOOLS
    echo
    echo "PATHTOOLS is set to $PATHTOOLS"
else
    echo "PATHTOOLS is already set to: $PATHTOOLS"
fi
export PATHTOOLS


if [ -z "$VMDPATH" ]; then
    # VMDPATH is not defined, so check if 'vmd' exists in the system
    if VMDPATH=$(command -v vmd); then
        echo "VMD found at: $VMDPATH"
        echo "VMDPATH is set to $VMDPATH"
    else
        echo "VMD not found in system"
#        VMDPATH=""
        exit
    fi
else
    echo "VMDPATH is already set to: $VMDPATH"
fi
export VMDPATH


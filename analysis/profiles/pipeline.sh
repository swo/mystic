#!/bin/bash

if [ ! -e "history.mat" ]; then
    echo history.mat is not here, make sure you copy it from interactive or wherever
    exit
fi

# make sure the plots folder exists
if [ ! -d "plots" ]; then
    echo Directory \"plots\" not found, making it
    mkdir plots
    mkdir plots/end
    mkdir plots/time
fi

matlab -nosplash -nodesktop -r 'make_plots; exit'

#!/bin/bash

if [ ! -e "history.mat" ]; then
    echo history.mat is not here, make sure you copy it from interactive or wherever
    exit
fi

matlab -nosplash -nodesktop -r 'get_divergences; exit'

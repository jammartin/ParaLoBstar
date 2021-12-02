#!/usr/bin/env bash

RUN_PATH="../../validation/plN4096th1_0"
RUNS="plN4096Se plN4096lb1np8Le plN4096lb1np8Hi"

for run in $RUNS; do
    ./main.py -d "$RUN_PATH/$run" # mass and energy
    ./main.py -d "$RUN_PATH/$run" -L # angular momentum
    ./main.py -d "$RUN_PATH/$run" -Q -l # mass quantiles
    mv output/energy_mass.png "output/${run}_energy_mass.png"
    mv output/angular_momentum.png "output/${run}_angular_momentum.png"
    mv output/mass_quantiles.png "output/${run}_mass_quantiles.png"
done

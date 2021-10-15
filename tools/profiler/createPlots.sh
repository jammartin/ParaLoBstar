#!/usr/bin/env bash

H5FILES="data/N1e6np8lb100CSGdxy5.h5 data/N1e6np8lb100CSGdxy3.h5 data/N1e6np8lb10CSGdxy3.h5"
H5PATHS_AMOUNT="general/numberOfParticles loadBalancing/sendParticles/receiveLength compF_BHpar/receiveLength"
H5PATHS_TIME="loadBalancing/totalTime forceComputation/totalTime updatePosVel/totalTime"

for file in $H5FILES; do
    for path_amount in $H5PATHS_AMOUNT; do
	./main.py --h5file "$file" --h5path "$path_amount"
    done
    for path_time in $H5PATHS_TIME; do
	./main.py --h5file "$file" --h5path "$path_time" -t
    done
done

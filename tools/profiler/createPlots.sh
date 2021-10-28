#!/usr/bin/env bash

H5FILES="data/profiling.h5"
H5PATHS_AMOUNT="general/numberOfParticles loadBalancing/receiveLength force/receiveLength"
H5PATHS_TIME="time/loadBalancing time/compForce"

for file in $H5FILES; do
    for path_amount in $H5PATHS_AMOUNT; do
	./main.py --h5file "$file" --h5path "$path_amount" -a
    done
    for path_time in $H5PATHS_TIME; do
	./main.py --h5file "$file" --h5path "$path_time" -t
    done
done

#!/usr/bin/env python3

import argparse
import sys
import numpy as np
import h5py


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Overview of execution time of paralobstar")
    parser.add_argument("--file", "-f",  metavar="str", type=str, help="path to the h5 profiling data file", required=True)
    parser.add_argument("--serial", "-s", action="store_true", help="analyzing output of the serial mode")
    parser.add_argument("--loadBalancing", "-l", metavar="int", type=int, help="load balancing interval (required if -s is not set)")
    parser.add_argument("--onlyTotal", "-o", action="store_true", help="only print total execution time per step in milliseconds")

    args = parser.parse_args()

    if not args.serial and not args.loadBalancing:
        print("Use -s for serial mode or provide load balancing interval! - Aborting.")
        sys.exit("Invalid combination of command line arguments. Either -s or -l required.")
    
    performance = h5py.File(args.file, 'r')

    # calculating averaged times per step
    pos = np.mean(np.amax(performance["/time/compPosition"], axis=1)) * 1e3
    move = np.mean(np.amax(performance["/time/moveParticles"], axis=1)) * 1e3
    if not args.serial:
        lb = np.mean(np.amax(performance["/time/loadBalancing"], axis=1)) * 1e3
        cc = np.mean(np.amax(performance["/time/commonCoarsePseudoParticles"], axis=1)) * 1e3
        fExchng = np.mean(np.amax(performance["/time/forceParticleExchange"], axis=1)) * 1e3
    pseudo = np.mean(np.amax(performance["/time/compPseudoParticles"], axis=1)) * 1e3
    force = np.mean(np.amax(performance["/time/compForce"], axis=1)) * 1e3
    vel = np.mean(np.amax(performance["/time/compVelocity"], axis=1)) * 1e3

    total = pos + move + pseudo + force + vel
    if not args.serial:
        total = total + lb/np.float64(args.loadBalancing)

    if args.onlyTotal:
        print(f'{total:.5f}')
    else:        
        print(f'{"OPERATION":30}| {"TIME [ms]":>14} | {"SHARE":>14} |')
        print("------------------------------|----------------|----------------|")
        print(f'{"compPosition":30}| {pos:14.3f} | {pos/total:14.2%} |')
        print(f'{"moveParticles":30}| {move:14.3f} | {move/total:14.2%} |')
        if not args.serial:
            print(f'{"loadBalancing":30}| {lb:14.3f} | {lb/args.loadBalancing/total:14.2%} | NOTE: is not done every time step, this is respected only in share')
            print(f'{"commonCoarsePseudoParticles":30}| {cc:14.3f} | {cc/total:14.2%} | NOTE: included in compPseudoParticles')
            print(f'{"forceParticleExchange":30}| {fExchng:14.3f} | {fExchng/total:14.2%} | NOTE: included in compForce')
        print(f'{"compPseudoParticles":30}| {pseudo:14.3f} | {pseudo/total:14.2%} |')
        print(f'{"compForce":30}| {force:14.3f} | {force/total:14.2%} |')
        print(f'{"compVelocity":30}| {vel:14.3f} | {vel/total:14.2%} |')
        print("------------------------------|----------------|----------------|")
        print(f'{"TOTAL":30}| {total:14.3f} |')


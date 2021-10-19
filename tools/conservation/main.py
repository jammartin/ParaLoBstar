#!/usr/bin/env python3

import argparse
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Plot conservation of energy and angular momentum for paralobstar's output.")
    parser.add_argument("--data", "-d", metavar="str", type=str, help="output data directory of paralobstar",
                        nargs="?", default="../../paralobstar/output")
    parser.add_argument("--angular_momentum", "-L", action="store_true", help="plot angular momentum (defaul: energy and mass)")

    args = parser.parse_args()

    time = []
    energy = []
    mass = []
    angular_momentum = []
    
    for h5file in sorted(glob.glob(os.path.join(args.data, "*.h5")), key=os.path.basename):
        print("Processing ", h5file, " ...")
        data = h5py.File(h5file, 'r')
        time.append(data["t"][()])
        energy.append(data["E_tot"][()])
        mass.append(np.sum(data["m"][:]))
        angular_momentum.append(np.array(data["L_tot"][:]))
        
        print("... done.")
                    
    plt.style.use("dark_background")
    fig, ax1 = plt.subplots(figsize=(12, 9), dpi=200)
    fig.patch.set_facecolor("black")
    ax1.set_xlabel("Time")
    
    if args.angular_momentum:
        ax1.set_title("Angular momentum")

        angMom = np.array(angular_momentum)
        
        ax1.plot(time, angMom[:,0], label="L_x")
        ax1.plot(time, angMom[:,1], label="L_y")
        ax1.plot(time, angMom[:,2], label="L_z")

        plt.legend(loc="best")
        
        fig.tight_layout()
        plt.savefig("output/angular_momentum.png")
        
    else:
        ax1.set_title("Total energy and mass")
        ax1.set_ylabel("Energy")
    
        ax1.plot(time, energy)

        ax2 = ax1.twinx()
        ax2.plot(time, mass)
        ax2.set_ylabel("Mass")

        fig.tight_layout()
        plt.savefig("output/energy_mass.png")

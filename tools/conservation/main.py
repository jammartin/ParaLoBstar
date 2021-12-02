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
    parser.add_argument("--mass_quantiles", "-Q", action="store_true", help="plot 10, 50 and 90 percent mass quantiles (default: energy and mass)")
    parser.add_argument("--light_theme", "-l", action="store_true", help="use a light theme for the generated plots")

    args = parser.parse_args()

    time = []
    energy = []
    mass = []
    angular_momentum = []
    mass_quantiles = []
    
    
    for h5file in sorted(glob.glob(os.path.join(args.data, "*.h5")), key=os.path.basename):
        print("Processing ", h5file, " ...")
        data = h5py.File(h5file, 'r')
        time.append(data["t"][()])

        if args.angular_momentum:
            print("... reading angular momentum ...")
            angular_momentum.append(np.array(data["L_tot"][:]))

        elif args.mass_quantiles:
            print("... computing mass quantiles ...")
            vecs2com = data["x"][:] - data["COM"][:]
            radii = np.linalg.norm(vecs2com, axis=1)
            radii.sort()
            numParticles = len(data["m"])
            print("NOTE: Only works for equal mass particle distributions!")
            mass_quantiles.append(np.array([
                radii[int(np.ceil(.1 * numParticles))],
                radii[int(np.ceil(.5 * numParticles))],
                radii[int(np.ceil(.9 * numParticles))]]))
        else:
            print("... computing mass and reading energy ...")
            energy.append(data["E_tot"][()])
            mass.append(np.sum(data["m"][:]))
        
        print("... done.")

        
    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{siunitx}')
    if not args.light_theme: plt.style.use("dark_background")
    fig, ax1 = plt.subplots(figsize=(8, 6), dpi=200)
    if not args.light_theme: fig.patch.set_facecolor("black")
    ax1.set_xlabel("Time $t$")
    
    if args.angular_momentum:
        ax1.set_title("Angular momentum along the three coordinate axes")

        angMom = np.array(angular_momentum)
        
        ax1.plot(time, angMom[:, 0], label="$L_x$")
        ax1.plot(time, angMom[:, 1], label="$L_y$")
        ax1.plot(time, angMom[:, 2], label="$L_z$")
        ax1.set_ylabel("Total angular momentum $L_c$ along axis $c$")
        
        plt.legend(loc="best")
        
        fig.tight_layout()
        plt.savefig("output/angular_momentum.png")

    elif args.mass_quantiles:
        ax1.set_title("Radii containing percentage of total mass")

        quantiles = np.array(mass_quantiles)

        ax1.plot(time, quantiles[:, 0], 'b', label="\SI{10}{\percent} mass quantile")
        ax1.plot(time, quantiles[:, 1], 'r', label="\SI{50}{\percent} mass quantile")
        if args.light_theme:
            ax1.plot(time, quantiles[:, 2], 'k', label="\SI{90}{\percent} mass quantile")
        else:
            ax1.plot(time, quantiles[:, 2], 'y', label="\SI{90}{\percent} mass quantile")
        
        ax1.set_ylabel("Radius $r$ containing a quantile of the total mass $M$")
        plt.legend(loc="best")
        plt.grid()
        fig.tight_layout()
        plt.savefig("output/mass_quantiles.png")


    else:
        ax1.set_title("Total energy and mass")
        ax1.set_ylabel("Total energy r$E_\text{tot}$")
    
        ax1.plot(time, energy, "r-", label=r"$E_\text{tot}$")

        ax2 = ax1.twinx()
        ax2.plot(time, mass, "b-", label="$M$")
        ax2.set_ylabel("Total mass $M$")

        fig.tight_layout()
        fig.legend()
        plt.savefig("output/energy_mass.png")

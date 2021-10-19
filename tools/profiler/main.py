#!/usr/bin/env python3

import argparse
import re
import numpy as np
#import matplotlib as mpl
import matplotlib.pyplot as plt
import h5py

SHARE_COLORS = ["#00ff00", "#0000ff", "#ff007f", "#ffff00",
                "#ff00ff", "#ff0000", "#00ffff", "#ff7f00"]
LINE_COLOR = "#5f5f5f"


def plotDataPerProcAsPercentage(h5file, path, axis, ylabel="Share per process"):
    data = np.array(h5file[path])
    data[data < 0] = 0  # removing all -1 values

    percDataBuff = np.transpose(data)
    # container for percentage data
    percData = np.zeros((data.shape[1], data.shape[0]), dtype=np.float64)

    for iProc in range(data.shape[1]):
        # Runtime Warning occurs in zeroth timestep when no particles are sent
        percData[iProc] = np.divide(percDataBuff[iProc, :].astype(np.float64), np.sum(data, 1).astype(np.float64))

    axis.set_ylabel(ylabel)
    axis.stackplot(np.arange(data.shape[0]), *percData, colors=SHARE_COLORS)


def plotTotalPerProc(h5file, path, axis, ylabel="Total per process"):
    data = np.array(h5file[path])
    data[data < 0] = 0  # removing all -1 values
    totalData = np.transpose(data)
    axis.set_ylabel(ylabel)
    axis.stackplot(np.arange(data.shape[0]), *totalData, colors=SHARE_COLORS)


def plotSumOverProcs(h5file, path, axis, ylabel="Sum over all processes"):
    data = h5file[path]
    axis.plot(np.arange(data.shape[0]), np.sum(data, 1), color=LINE_COLOR)
    axis.set_ylabel(ylabel)


def plotMaxOfProcs(h5file, path, axis, ylabel="Maximum of all processes"):
    data = h5file[path]
    axis.plot(np.arange(data.shape[0]), np.max(data, 1), color=LINE_COLOR)
    axis.set_ylabel(ylabel)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Plot profifling data for paralobstar.")
    parser.add_argument("--h5path", metavar="str", type=str, help="h5 file path to the data to plot",
                        nargs="?", default="general/numberOfParticles")
    parser.add_argument("--time", "-t", action="store_true", help="plotting time instead of amounts")
    parser.add_argument("--h5file", metavar="str", type=str, help="input h5 file path",
                        nargs="?", default="data/performance.h5")
    parser.add_argument("--absolute", "-a", action="store_true", help="plot absolute stack without total")

    args = parser.parse_args()

    print("Plotting data from '" + args.h5file + "'@/" + args.h5path + "...")

    performance = h5py.File(args.h5file, 'r')

    # Initialize plot style
    plt.style.use('dark_background')


    fig, ax1 = plt.subplots(figsize=(12, 9), dpi=200)
    fig.patch.set_facecolor("black")
    ax1.set_title("/" + args.h5path)
    ax1.set_xlabel("Timestep")

    if args.absolute:
        plotTotalPerProc(performance, args.h5path, ax1)

    else:
        plotDataPerProcAsPercentage(performance, args.h5path, ax1)

        ax2 = ax1.twinx()

        if args.time:
            plotMaxOfProcs(performance, args.h5path, ax2)
        else:
            plotSumOverProcs(performance, args.h5path, ax2)

    fig.tight_layout()
    plt.savefig(re.sub(r"data/(.*)\.h5", r"output/\g<1>_" + args.h5path.replace("/", "-") + ".png", args.h5file))

    print("... done.")

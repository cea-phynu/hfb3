#!/usr/bin/env python3
# This script can be used to gather some observables from a PES into a .csv file.

# ==============================================================================
# ==============================================================================
# ==============================================================================


import pickle
import numpy as np
import logging

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s.%(msecs)03d [%(levelname)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

# ==============================================================================
# ==============================================================================
# ==============================================================================


def genCsv(path, beta2min, beta2max, beta2step):
    pickleFile = path + "/checkpoint.pickle"
    try:
        data = pickle.load(open(pickleFile, "rb"))
        points, history, plotId = data

        nbPoints = 0
        for p in points:
            if p["md5sum"]:
                nbPoints += 1

        logging.info(f"found checkpoint ({len(history)} hist. {nbPoints} points)")

    except Exception:
        return

    beta2   = [round(float(beta2), 2) for beta2 in np.arange(beta2min, beta2max + .001, beta2step)]

    # find minimum energy
    eneMin = 1e99
    for id in range(len(points)):
        if points[id]['totalEnergy'] < eneMin and abs(beta2[id]) < 0.55:
            eneMin = points[id]['totalEnergy']

    # generate .csv file
    path = path + "_normalized.csv"
    logging.info(f"saving in {path}")

    with open(path, "w+") as fp:
        fp.write("beta2,Ehfb,Ehfb normalized\n")
        for id in range(len(points)):
            if points[id]["md5sum"]:
                fp.write(f"{beta2[id]},{points[id]['totalEnergy']},{points[id]['totalEnergy'] - eneMin}\n")

# ==============================================================================
# ==============================================================================
# ==============================================================================


if __name__ == "__main__":

    import sys
    if len(sys.argv) < 4:
        print(f"Usage: {sys.argv[0]} path beta2min beta2max beta2step")
        sys.exit(0)

    beta2min = float(sys.argv[2])
    beta2max = float(sys.argv[3])
    beta2step = float(sys.argv[4])
    genCsv(sys.argv[1], beta2min, beta2max, beta2step)

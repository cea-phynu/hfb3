#!/usr/bin/env python3
# This script can be used to gather some observables from a PES into a .csv file.

# ==============================================================================
# ==============================================================================
# ==============================================================================


import sys
import logging
import hfb3
import numpy as np

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s.%(msecs)03d [%(levelname)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

# ==============================================================================
# ==============================================================================
# ==============================================================================


def process(fileName):

    logging.info(f"{fileName}")

    d = hfb3.DataTree.getDefault() + hfb3.DataTree(fileName)
    s = hfb3.State(d)
    multipoleOperators = hfb3.MultipoleOperators(s)
    s.calcUVFromRhoKappa(d)
    # a = hfb3.Action(d)
    # a.calcObservables()

    coords = np.array((2,), dtype=np.int32  , order='F')
    s.calcInertia(coords)
    result = {}
    result["beta2"]       = multipoleOperators.beta[hfb3.TOTAL]
    result["q20"]         = multipoleOperators.qlm[2, hfb3.TOTAL]
    result["q30"]         = multipoleOperators.qlm[3, hfb3.TOTAL]
    result["energy"]      = s.totalEnergy
    result["converged"]   = 1 if s.converged else 0
    result["zpeATDHF00"]  = s.zpeATDHF[0, 0]
    result["massATDHF00"] = s.massATDHF[0, 0]
    result["zpeGCM"]      = s.zpeGCM[0, 0]
    result["massGCM"]     = s.massGCM[0, 0]
    # result.append(fileName)
    return result

# ==============================================================================
# ==============================================================================
# ==============================================================================


if __name__ == "__main__":

    from multiprocessing import Pool

    path = "."

    if len(sys.argv) > 1:
        files = sys.argv[1:]

    hfb3.cvar.msgToOut = (hfb3.MSG_ERROR,)

    p = Pool(11)
    results = p.map(process, files)

    final = {}

    columnNames = []
    for r in results:
        beta2 = round(r["beta2"], 3)
        if r["converged"] == 0:
            continue
        if beta2 not in final:
            final[beta2] = r
        if r["energy"] < final[beta2]["energy"]:
            final[beta2] = r
        columnNames = [c for c in r.keys()]

    sortedResults = sorted(final.values(), key=lambda v: v["beta2"])

    print(",".join(columnNames))

    for r in sortedResults:
        ifirst = True
        for f in columnNames:
            if ifirst:
                ifirst = False
            else:
                print(",", end='')
            print(f"{r[f]}", end='')
        print("")

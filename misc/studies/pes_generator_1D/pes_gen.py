#!/usr/bin/env python3
# Study of the influence of the basis size on the HFB energy for a set of
# nuclei.

# ==============================================================================
# ==============================================================================
# ==============================================================================


import hfb3
import os
import pickle
import numpy as np
import multiprocessing
import hashlib
import time
import logging
import psutil
import gzip
import shutil

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s.%(msecs)03d [%(levelname)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

# ==============================================================================
# ==============================================================================
# ==============================================================================


def setCPUAdffinity():
    # fix CPU affinity (messed up on some calculators)
    os.unsetenv('KMP_AFFINITY')
    nbCores = os.cpu_count()
    affinityMask = set(range(nbCores))
    os.sched_setaffinity(0, affinityMask)

# ==============================================================================
# ==============================================================================
# ==============================================================================


def worker(workDir, resultsQueue, dataTree, beta2, fromMD5):

    class MyCallBack(hfb3.Callback):

        def __init__(self, filename):
            self.filename = filename
            hfb3.Callback.__init__(self)

        def run(self, str):
            with open(self.filename, "a+") as fp:
                fp.write(str + "\n")
                # print(str, flush = True)

    if fromMD5:
        initialState = workDir + "/" + fromMD5 + ".msg.gz"
        dataTree.merge(hfb3.DataTree(initialState))
        inputFile = "%5.3f_%s" % (beta2, fromMD5)
    else:
        inputFile = "%5.3f_scratch" % (beta2,)

    outputFile = workDir + "/" + inputFile + ".out"
    hfb3.cvar.logger.setCallback(MyCallBack(outputFile).__disown__())

    dataTree.setD("constraints/beta2t", beta2)
    dataTree.setB("action/saveResultFiles", False)
    action = hfb3.Action(dataTree)
    action.calcHFB()

    content = hfb3.dataTreeToBytes(action.state.getDataTree())
    md5sum = hashlib.md5(content).hexdigest()

    resultFile = workDir + "/" + md5sum + ".msg.gz"
    with gzip.open(resultFile, "wb+") as fp:
        fp.write(content)
        hfb3.Tools.info("result saved in " + resultFile)

    resultsQueue.put((beta2, action.state.totalEnergy, action.state.converged, md5sum, fromMD5))

# ==============================================================================
# ==============================================================================
# ==============================================================================


class Pes:

    # ==========================================================================

    def __init__(self, beta2min, beta2max, beta2step, hfb3File, label=""):

        setCPUAdffinity()

        self.name = hfb3File + label
        self.workDir = hfb3File + label + ".temp"
        self.hfb3File = hfb3File

        try:
            os.mkdir(self.workDir)
        except Exception:
            pass

        self.pickleFile = self.workDir + "/checkpoint.pickle"

        defaultDataTree = hfb3.DataTree.getDefault()
        self.dataTree = defaultDataTree + hfb3.DataTree(hfb3File)

        self.resultsQueue = multiprocessing.Queue()

        self.beta2   = [round(float(beta2), 2) for beta2 in np.arange(beta2min, beta2max + .001, beta2step)]
        self.points  = [dict(totalEnergy=1e99, md5sum=None) for b in self.beta2]
        self.history = dict()

        self.nbRunning = 0
        self.plotId = 0

        self.restore()

    # ==========================================================================

    def backup(self):

        nbPoints = 0
        for p in self.points:
            if p["md5sum"]:
                nbPoints += 1

        historyCleaned = {}
        for k, v in self.history.items():
            if v is not None :
                historyCleaned[k] = v

        data = [self.points, historyCleaned, self.plotId]
        pickle.dump(data, open(self.pickleFile, "wb+"))
        logging.info(f"{self.name}: saved checkpoint ({len(historyCleaned)} hist. {nbPoints} points) {self.nbRunning + 1}/{self.maxRunning}")

    # ==========================================================================

    def restore(self):
        try:
            data = pickle.load(open(self.pickleFile, "rb"))
            self.points, self.history, self.plotId = data

            nbPoints = 0
            for p in self.points:
                if p["md5sum"]:
                    nbPoints += 1

            logging.info(f"{self.name}: found checkpoint ({len(self.history)} hist. {nbPoints} points)")

            toDelete = []
            for k, v in self.history.items():
                if v is None:
                    toDelete.add(k)

            for k in toDelete:
                del self.history[k]

        except Exception:
            pass

    # ==========================================================================

    def launch(self, fromId, id):

        # logging.debug(f"{self.name}: trying to launch {fromId} {id}")
        if self.nbRunning >= self.maxRunning:
            return False

        if psutil.virtual_memory().percent > 90.0 and self.nbRunning > 0:
            return False

        fromMD5 = None
        if fromId is not None:
            fromMD5 = self.points[fromId]["md5sum"]
            if fromMD5 is None:
                return False

        key = (id, fromMD5)
        if key not in self.history.keys():
            # logging.debug(f"{self.name}: key not in history: {key}")
            multiprocessing.Process(target=worker, args=(self.workDir, self.resultsQueue, self.dataTree, self.beta2[id], fromMD5)).start()
            self.history[key] = None
            self.nbRunning += 1
            if fromId:
                logging.debug(f"{self.name}: launch {id}({self.beta2[id]}) from {fromId}({self.beta2[fromId]}) [{fromMD5}]")
            else:
                logging.debug(f"{self.name}: launch {id}({self.beta2[id]}) from scratch")
            return True

        return False

    # ==========================================================================

    def calc(self, nbProcs=4):

        self.maxRunning = nbProcs

        logging.info(f"{self.name}: calculating PES for {self.hfb3File} using {self.maxRunning} processes")

        firstTime = True

        while True:

            newResults = False

            # receive results
            while not self.resultsQueue.empty():
                result = self.resultsQueue.get()
                self.nbRunning -= 1

                beta2, totalEnergy, converged, md5sum, fromMD5 = result
                id = self.beta2.index(beta2)
                # logging.debug(f"{self.name}: result {beta2} {totalEnergy} ({md5sum})")

                key = (id, fromMD5)
                self.history[key] = md5sum

                oldEnergy = self.points[id]["totalEnergy"]
                if oldEnergy - 0.001 > totalEnergy:
                    self.points[id]["md5sum"] = md5sum
                    self.points[id]["totalEnergy"] = totalEnergy
                    logging.debug(f"{self.name}: update {id}({beta2}): {oldEnergy} -> {totalEnergy} [{md5sum}]")
                    self.printCsv(f"result_{self.plotId:06d}.csv")
                    shutil.copyfile(self.workDir + "/" + f"result_{self.plotId:06d}.csv", self.workDir + "/" + "result.csv")
                    self.plotId += 1

                newResults = True

            if newResults:
                # backup the PES
                self.backup()

            if newResults or firstTime:
                firstTime = False

                launchedAny = False
                launchedAny |= self.propag(0)  # initial mesh
                launchedAny |= self.propag(8)
                launchedAny |= self.propag(-8)
                launchedAny |= self.propag(4)
                launchedAny |= self.propag(-4)
                launchedAny |= self.propag(2)
                launchedAny |= self.propag(-2)
                launchedAny |= self.propag(1)
                launchedAny |= self.propag(-1)

            # go easy on the loop
            time.sleep(0.1)

            # continue ?
            if self.nbRunning == 0:
                break

    # ==========================================================================

    def propag(self, step):
        launchedAny = False
        if step > 0:
            for id in range(step, len(self.beta2), step):
                launchedAny |= self.launch(id - step, id)
        elif step < 0:
            for id in range(0, len(self.beta2) + step, -step):
                launchedAny |= self.launch(id - step, id)
        else:
            for id in range(0, len(self.beta2), 8):
                launchedAny |= self.launch(None, id)

        return launchedAny

    # ==========================================================================

    def printCsv(self, fileName="result.csv"):

        # find minimum energy
        eneMin = 1e99
        for id in range(len(self.points)):
            if self.points[id]['totalEnergy'] < eneMin and abs(self.beta2[id]) < 0.55:
                eneMin = self.points[id]['totalEnergy']

        # generate .csv file
        with open(self.workDir + "/" + fileName, "w+") as fp:
            fp.write("beta2,Ehfb,Ehfb normalized\n")
            for id in range(len(self.points)):
                if self.points[id]["md5sum"]:
                    fp.write(f"{self.beta2[id]},{self.points[id]['totalEnergy']},{self.points[id]['totalEnergy'] - eneMin}\n")

    # ==========================================================================

    def __str__(self):
        s = ""
        s += f'nb points  : {len(self.points)}'
        s += f'\nmin  beta20: {self.beta2[0]}'
        s += f'\nmax  beta20: {self.beta2[-1]}'
        s += f'\nstep beta20: {self.beta2[1] - self.beta2[0]}'
        return s

# ==============================================================================
# ==============================================================================
# ==============================================================================


if __name__ == "__main__":

    import sys
    if len(sys.argv) < 9:
        print(f"Usage: {sys.argv[0]} file.hfb3 beta2min beta2max beta2step centers nOscil g_q maxProcs")
        sys.exit(0)

    beta2min = float(sys.argv[2])
    beta2max = float(sys.argv[3])
    beta2step = float(sys.argv[4])
    centers = int(sys.argv[5])
    nOscil = int(sys.argv[6])
    g_q = float(sys.argv[7])
    maxProcs = int(sys.argv[8])

    label = f"_{centers}x{nOscil:02}_Q{g_q:.2f}"

    pes = Pes(beta2min=beta2min, beta2max=beta2max, beta2step=beta2step, hfb3File=sys.argv[1], label=label)
    pes.dataTree.setI("basis/nOscil", nOscil)
    pes.dataTree.setD("basis/g_q", g_q)
    if centers == 1:
        pes.dataTree.setD("basis/d_0", 0.0)
    pes.calc(maxProcs)

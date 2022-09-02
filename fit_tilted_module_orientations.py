#!/bin/env python

import ROOT as r
import math
import sys

import Module as m
import tqdm

f = open("hits.txt")
lines = f.readlines()
hits = {}
moduleType = {}

# Parse and store hits
for line in tqdm.tqdm(lines):
    if "detId" not in line:
        continue
    ls = line.split()
    detid = ls[7]
    if detid not in moduleType:
        moduleType[detid] = ls[9]
    hit = (float(ls[1]), float(ls[3]), float(ls[5])) # NOTE in txt file we have z, y, x coordinates 
    if detid not in hits:
        hits[detid] = []
    hits[detid].append(hit)

output = open("data/tilted_orientation.txt", "w")

# Computing two groups of hits
output.write("# detid drdz xy-slope\n")
for detid in tqdm.tqdm(hits):

    mod = m.Module(int(detid), int(moduleType[detid]))
    istilt = (mod.side() == 1 or mod.side() == 2) and mod.subdet() == 5
    if not istilt:
        continue
    isstrip = mod.moduleLayerType() == 1
    if not isstrip:
        continue

    # Number of events
    n = len(hits[detid])

    # There are two groups of hits

    azs = []
    for ii, hit in enumerate(hits[detid]):
        azs.append(abs(hit[2]))
    azs = list(set(azs))
    azs.sort()

    # Loop again and group them into two groups
    hl = [] # low-z hits
    hh = [] # high-z hits
    xls = []
    xhs = []
    for ii, hit in enumerate(hits[detid]):
        az = abs(hit[2])
        # put them into two buckets
        if az == azs[0]:
            hl.append(hit)
            xls.append(hit[0])
        else:
            hh.append(hit)
            xhs.append(hit[0])

    # Create two TGraph's for fit
    gl = r.TGraph(len(hl))
    gh = r.TGraph(len(hh))

    for ii, hit in enumerate(hl):
        gl.SetPoint(ii, hit[0], hit[1])

    # if lying 90 degrees in x-y plane the fit will fail with infinite slope
    # so take care of it as a special case

    if len(list(set(xls))) != 1:

        rl = gl.Fit("pol1", "q") # Result of low hits fit
        yl = r.gROOT.FindObject("pol1").GetParameter(0)
        sl = r.gROOT.FindObject("pol1").GetParameter(1)

        for ii, hit in enumerate(hh):
            gh.SetPoint(ii, hit[0], hit[1])

        rh = gh.Fit("pol1", "q") # Result of high hits fit
        yh = r.gROOT.FindObject("pol1").GetParameter(0)
        sh = r.gROOT.FindObject("pol1").GetParameter(1)

        if abs(sl - sh) > 0.005:
            print("ERROR")

        if abs(yh-yl)/math.sqrt(sl**2+1) == 0:
            print("")
            for h in hl:
                print(h)
            print("")
            for h in hh:
                print(h)
            rl = gl.Fit("pol0", "q") # Result of low hits fit
            yl = r.gROOT.FindObject("pol0").GetParameter(0)
            print(yl)
            sys.exit()

        output.write("{} {} {}\n".format(detid, abs(yh-yl)/math.sqrt(sl**2+1)/abs(azs[0]-azs[1]), sl)) #, abs(yh-yl)/math.sqrt(sl**2+1), abs(azs[0] - azs[1])

    else:

        output.write("{} {} {}\n".format(detid, abs(list(set(xls))[0]-list(set(xhs))[0])/ abs(azs[0] - azs[1]), 123456789)) #, abs(list(set(xls))[0]-list(set(xhs))[0]), abs(azs[0] - azs[1])


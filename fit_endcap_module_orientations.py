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

output = open("data/endcap_orientation.txt", "w")

# Computing two groups of hits
output.write("# detid average_r2s y_intercept_low_hits slope_low_hits y_intercept_high_hits slope_high_hits\n")
for detid in tqdm.tqdm(hits):

    mod = m.Module(int(detid), int(moduleType[detid]))
    isendcap = mod.subdet() == 4
    if not isendcap:
        continue
    isstrip = mod.moduleLayerType() == 1
    if not isstrip:
        continue

    # Number of events
    n = len(hits[detid])

    # There are two groups of hits
    # The average value of r^2's will be used to divide the hits into two groups

    # Trying to group into two hits
    r2s = []

    # compute r^2 = x^2 + y^2
    sumr2s = 0
    for ii, hit in enumerate(hits[detid]):
        r2 = hit[0]**2 + hit[1]**2
        r2s.append(r2)
        sumr2s += r2

    # The average value
    avgr2s = sumr2s / n

    # Loop again and group them into two groups
    hl = [] # low hits
    hh = [] # high hits
    for ii, hit in enumerate(hits[detid]):
        r2 = hit[0]**2 + hit[1]**2
        # put them into two buckets
        if r2 < avgr2s:
            hl.append(hit)
        else:
            hh.append(hit)

    # Create two TGraph's for fit
    gl = r.TGraph(len(hl))
    gh = r.TGraph(len(hh))

    for ii, hit in enumerate(hl):
        gl.SetPoint(ii, hit[0], hit[1])

    if len(hl) > 1:
        rl = gl.Fit("pol1", "q") # Result of low hits fit
        yl = r.gROOT.FindObject("pol1").GetParameter(0)
        sl = r.gROOT.FindObject("pol1").GetParameter(1)
    else:
        rl = -999
        yl = -999
        sl = -999

    for ii, hit in enumerate(hh):
        gh.SetPoint(ii, hit[0], hit[1])

    if len(hh) > 1:
        rh = gh.Fit("pol1", "q") # Result of high hits fit
        yh = r.gROOT.FindObject("pol1").GetParameter(0)
        sh = r.gROOT.FindObject("pol1").GetParameter(1)
    else:
        rh = -999
        yh = -999
        sh = -999

    # if abs(sl - sh) > 0.5 and (sl != -999 and sh != -999):
    #     print("ERROR", sl, sh)
    #     print(hl)
    #     print(hh)

    #     for i in hl:
    #         print(i)

    #     print("")
    #     for i in hh:
    #         print(i)

    if sl == -999 and sh != -999:
        sl = sh
        yl = yh
        rl = rh

    if sl != -999 and sh == -999:
        sh = sl
        yh = yl
        rh = rl

    # if sl == -999:
    #     print(detid, avgr2s, yl, sl, yh, sh)

    output.write("{} {} {} {} {} {}\n".format(detid, avgr2s, yl, sl, yh, sh))


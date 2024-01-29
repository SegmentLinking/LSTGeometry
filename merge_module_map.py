#!/bin/env python

import sys
import os
from tqdm import tqdm

connections = {}

# Take file `f` and `g` and combine into `h`

f = open("output/module_connection_tracing_straight.txt")
g = open("output/module_connection_tracing_curved.txt")
h = open("output/module_connection_tracing_merged.txt", "w")

lines = f.readlines()

for line in tqdm(lines, desc="Looping over 1st connection map"):
    ls = line.split()
    ref_detid = int(ls[0])
    nconn = int(ls[1])
    target_detids = [ int(x) for x in ls[2:] ]

    if ref_detid not in connections:
        connections[ref_detid] = []

    connections[ref_detid] += target_detids

lines = g.readlines()

for line in tqdm(lines, desc="Looping over 2nd connection map"):
    ls = line.split()
    ref_detid = int(ls[0])
    nconn = int(ls[1])
    target_detids = [ int(x) for x in ls[2:] ]

    if ref_detid not in connections:
        connections[ref_detid] = []

    connections[ref_detid] += target_detids

for reference in sorted(tqdm(connections.keys())):
    uniquelist = list(set(connections[reference]))
    nconn = len(uniquelist)
    targets = [ str(x) for x in uniquelist ]
    h.write("{} {} {}\n".format(reference, nconn, " ".join(targets)))

#!/bin/env python

import tqdm

f = open("./hits.txt")
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

output = open("data/centroid.txt", "w")

for detid in hits:

    nhits = len(hits[detid])
    sum_x = 0
    sum_y = 0
    sum_z = 0
    for hit in hits[detid]:
        sum_x += hit[0]
        sum_y += hit[1]
        sum_z += hit[2]
    avg_x = sum_x / nhits
    avg_y = sum_y / nhits
    avg_z = sum_z / nhits

    output.write("{},{},{},{},{}\n".format(detid, avg_x, avg_y, avg_z, moduleType[detid]))

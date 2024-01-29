# Goal is to compute something like the following

# # average values
# self.average_radii = []
# self.average_radii.append(24.726409077007705) # Layer 1 average radius
# self.average_radii.append(37.059873804403495) # Layer 2 average radius
# self.average_radii.append(52.17677700048082)  # Layer 3 average radius
# self.average_radii.append(68.61016946477243)  # Layer 4 average radius
# self.average_radii.append(85.91013998484999)  # Layer 5 average radius
# self.average_radii.append(110.71009476599565) # Layer 6 average radius
# self.average_zs = []
# self.average_zs.append(130.93374689440995) # Layer 1 average Z (endcap)
# self.average_zs.append(154.74990605590062) # Layer 2 average Z (endcap)
# self.average_zs.append(185.1167890070922)  # Layer 3 average Z (endcap)
# self.average_zs.append(221.39607712765957) # Layer 4 average Z (endcap)
# self.average_zs.append(264.76252304964544) # Layer 5 average Z (endcap)

import os
import sys
# Add the parent directory to sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import tqdm
import math
from Module import Module

f = open("hits.txt")
lines = f.readlines()
hits = {}

# Parse and store hits
for line in tqdm.tqdm(lines):
    if "detId" not in line:
        continue
    ls = line.split()
    detid = ls[7]
    hit = (float(ls[1]), float(ls[3]), float(ls[5])) # NOTE in txt file we have z, y, x coordinates 
    if detid not in hits:
        hits[detid] = []
    hits[detid].append(hit)

average_radius_output = open("../output/average_radius.txt", "w")
average_z_output = open("../output/average_z.txt", "w")

radii_layer_1_barrel = []
radii_layer_2_barrel = []
radii_layer_3_barrel = []
radii_layer_4_barrel = []
radii_layer_5_barrel = []
radii_layer_6_barrel = []

z_layer_1_endcap = []
z_layer_2_endcap = []
z_layer_3_endcap = []
z_layer_4_endcap = []
z_layer_5_endcap = []

for detid in hits:

    mod = Module(detid)

    isBarrel = mod.subdet() == 5
    layer = mod.layer()

    nhits = len(hits[detid])
    for hit in hits[detid]:
        x = hit[0]
        y = hit[1]
        z = hit[2]
        r = math.sqrt(x**2 + y**2)

        if isBarrel and layer == 1: radii_layer_1_barrel.append(r)
        if isBarrel and layer == 2: radii_layer_2_barrel.append(r)
        if isBarrel and layer == 3: radii_layer_3_barrel.append(r)
        if isBarrel and layer == 4: radii_layer_4_barrel.append(r)
        if isBarrel and layer == 5: radii_layer_5_barrel.append(r)
        if isBarrel and layer == 6: radii_layer_6_barrel.append(r)

        if not isBarrel and layer == 1: z_layer_1_endcap.append(abs(z))
        if not isBarrel and layer == 2: z_layer_2_endcap.append(abs(z))
        if not isBarrel and layer == 3: z_layer_3_endcap.append(abs(z))
        if not isBarrel and layer == 4: z_layer_4_endcap.append(abs(z))
        if not isBarrel and layer == 5: z_layer_5_endcap.append(abs(z))

avg_radius_layer_1_barrel = sum(radii_layer_1_barrel) / len(radii_layer_1_barrel)
avg_radius_layer_2_barrel = sum(radii_layer_2_barrel) / len(radii_layer_2_barrel)
avg_radius_layer_3_barrel = sum(radii_layer_3_barrel) / len(radii_layer_3_barrel)
avg_radius_layer_4_barrel = sum(radii_layer_4_barrel) / len(radii_layer_4_barrel)
avg_radius_layer_5_barrel = sum(radii_layer_5_barrel) / len(radii_layer_5_barrel)
avg_radius_layer_6_barrel = sum(radii_layer_6_barrel) / len(radii_layer_6_barrel)

avg_z_layer_1_endcap = sum(z_layer_1_endcap) / len(z_layer_1_endcap)
avg_z_layer_2_endcap = sum(z_layer_2_endcap) / len(z_layer_2_endcap)
avg_z_layer_3_endcap = sum(z_layer_3_endcap) / len(z_layer_3_endcap)
avg_z_layer_4_endcap = sum(z_layer_4_endcap) / len(z_layer_4_endcap)
avg_z_layer_5_endcap = sum(z_layer_5_endcap) / len(z_layer_5_endcap)

average_radius_output.write("{}\n".format(avg_radius_layer_1_barrel))
average_radius_output.write("{}\n".format(avg_radius_layer_2_barrel))
average_radius_output.write("{}\n".format(avg_radius_layer_3_barrel))
average_radius_output.write("{}\n".format(avg_radius_layer_4_barrel))
average_radius_output.write("{}\n".format(avg_radius_layer_5_barrel))
average_radius_output.write("{}\n".format(avg_radius_layer_6_barrel))

average_z_output.write("{}\n".format(avg_z_layer_1_endcap))
average_z_output.write("{}\n".format(avg_z_layer_2_endcap))
average_z_output.write("{}\n".format(avg_z_layer_3_endcap))
average_z_output.write("{}\n".format(avg_z_layer_4_endcap))
average_z_output.write("{}\n".format(avg_z_layer_5_endcap))

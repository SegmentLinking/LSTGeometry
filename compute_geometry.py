#!/bin/env python

# Script to compute the detector module rectangle plane (the 4 corners) from the centroid data

# Output format
# {
#     "411309061": [
#         [
#             -129.13503188311316,   # z coord
#             27.99223639171157,     # y coord
#             -0.7045785443457637    # x coord
#         ],
#         [
#             -129.13503188311316,
#             26.454739008300713,
#             9.176519663647351
#         ],
#         [
#             -129.13503188311316,
#             21.514189904304153,
#             8.407770971941925
#         ],
#         [
#             -129.13503188311316,
#             23.051687287715012,
#             -1.4733272360511913
#         ]
#     ],
#     ...
#     ...
#     ...
#     ...
#     ...
# }

from Centroid import Centroid
from TiltedOrientation import TiltedOrientation
import Module as m

import math
import json

centroid = Centroid("data/sensor_centroids.txt")
tiltedorientation = TiltedOrientation("data/tilted_orientation.txt")

module_four_corners_database = {}

for index, detid in enumerate(centroid.data):

    x, y, z, moduleTypeInfo = centroid.getCentroid(detid)

    module = m.Module(detid, moduleTypeInfo)

    # if module.side() != 1:
    #     continue
    # if module.layer() != 1 and module.layer() != 6:
    #     continue
    # if module.rod() != 12:
    #     continue
    # if module.subdet() == 5:
    #     continue
    # if module.layer() != 2:
    #     continue
    # if module.module() != 1:
    #     continue
    # if module.moduleType() == 1:
    #     continue

    r = math.sqrt(x**2 + y**2)
    phi = math.atan2(y, x)
    phi_ = math.atan2(x, -y)
    is2S = module.moduleType()
    philen = 5
    etalen = 5 if is2S else 2.5

    # Barrel
    if module.subdet() == 5:

        # Either drdz = None or a value
        drdz = None
        if module.subdet() == 5 and module.side() != 3: # If it is tilted
            try:
                drdz = tiltedorientation.getDrDz(detid)
            except:
                drdz = tiltedorientation.getDrDz(module.partnerDetId())

        # Compute various variables
        dz = etalen
        dr = drdz * dz if drdz else 0
        norm = math.sqrt(dz**2 + dr**2)
        dz /= norm
        dr /= norm
        dz *= etalen
        dr *= etalen
        dx = philen * math.cos(phi_)
        dy = philen * math.sin(phi_)
        dxr = dr * math.cos(phi)
        dyr = dr * math.sin(phi)

        # Change sign of tilt related displacement if other side of the hemisphere
        if module.side() == 2 and module.subdet() == 5:
            dxr = -dxr
            dyr = -dyr

    # Endcap
    else:

        # Compute various variables
        dz = 0
        dr = etalen
        norm = math.sqrt(dz**2 + dr**2)
        dz /= norm
        dr /= norm
        dz *= etalen
        dr *= etalen
        dx = philen * math.cos(phi_)
        dy = philen * math.sin(phi_)
        dxr = dr * math.cos(phi)
        dyr = dr * math.sin(phi)

    # Compute the four corners
    four_corner_coords = [
        [ z + dz , x + dx + dxr , y + dy + dyr ],
        [ z - dz , x + dx - dxr , y + dy - dyr ],
        [ z - dz , x - dx - dxr , y - dy - dyr ],
        [ z + dz , x - dx + dxr , y - dy + dyr ],
    ]

    # push to the list
    module_four_corners_database[detid] = four_corner_coords

g = open("data/sensor_corners.txt", "w")
g.write(json.dumps(module_four_corners_database, indent=4))

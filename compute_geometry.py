import sys
import json
import numpy as np
import Module as m
from Centroid import Centroid
from TiltedOrientation import TiltedOrientation

# Compute the detector module rectangle plane (the 4 corners) from the centroid data

# Output format
# {
#     "411309061": [
#         [
#             -129.13503188311316,
#             27.99223639171157,
#             -0.7045785443457637
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
#     ], ...
# }

def compute_four_corners_barrel(x, y, z, phi, phi_, etalen, philen, drdz, module_side):
    """
    Computes the four corner coordinates for a barrel module.

    Parameters:
    - x, y, z: Centroid coordinates
    - phi, phi_: Angles in the xy-plane
    - etalen: Length along the eta direction
    - philen: Length along the phi direction
    - drdz: Slope of the module
    - module_side: Side of the module (1, 2, or 3)

    Returns:
    - List of four corner coordinates
    """
    dz = etalen
    dr = drdz * dz if drdz else 0
    norm = np.sqrt(dz**2 + dr**2)
    dz /= norm
    dr /= norm
    dz *= etalen
    dr *= etalen
    dx = philen * np.cos(phi_)
    dy = philen * np.sin(phi_)
    dxr = dr * np.cos(phi)
    dyr = dr * np.sin(phi)

    if module_side == 2:
        dxr = -dxr
        dyr = -dyr

    return [
        [z + dz, x + dx + dxr, y + dy + dyr],
        [z - dz, x + dx - dxr, y + dy - dyr],
        [z - dz, x - dx - dxr, y - dy - dyr],
        [z + dz, x - dx + dxr, y - dy + dyr]
    ]

def compute_four_corners_endcap(x, y, z, phi, phi_, etalen, philen):
    """
    Computes the four corner coordinates for an endcap module.

    Parameters:
    - x, y, z: Centroid coordinates
    - phi, phi_: Angles in the xy-plane
    - etalen: Length along the eta direction
    - philen: Length along the phi direction

    Returns:
    - List of four corner coordinates
    """
    dz = 0
    dr = etalen
    norm = np.sqrt(dz**2 + dr**2)
    dz /= norm
    dr /= norm
    dz *= etalen
    dr *= etalen
    dx = philen * np.cos(phi_)
    dy = philen * np.sin(phi_)
    dxr = dr * np.cos(phi)
    dyr = dr * np.sin(phi)

    return [
        [z + dz, x + dx + dxr, y + dy + dyr],
        [z - dz, x + dx - dxr, y + dy - dyr],
        [z - dz, x - dx - dxr, y - dy - dyr],
        [z + dz, x - dx + dxr, y - dy + dyr]
    ]

def get_drdz(tiltedorientation, detid, partner_id):
    """
    Safely retrieves the drdz value for a given detid.
    If not found, tries the partner_id.
    """
    if detid in tiltedorientation.data:
        return tiltedorientation.getDrDz(detid)
    if partner_id in tiltedorientation.data:
        return tiltedorientation.getDrDz(partner_id)
    return None

if __name__ == "__main__":
    # Default file paths
    default_centroid_file = "data/centroid.txt"
    default_orientation_file = "data/tilted_orientation.txt"
    default_output_file = "data/geom2.txt"

    # Check for help flag
    if '-h' in sys.argv or '--help' in sys.argv:
        print("\nUsage: python compute_geometry.py [centroid_file] [orientation_file] [output_file]")
        print("\nOptions:")
        print(f"  centroid_file      Path to the input file containing centroid data. Default is {default_centroid_file}")
        print(f"  orientation_file   Path to the input file containing tilted orientation data. Default is {default_orientation_file}")
        print(f"  output_file        Path for the output file. Default is {default_output_file}\n")
        sys.exit()

    # Determine file paths based on arguments provided
    centroid_file = sys.argv[1] if len(sys.argv) > 1 else default_centroid_file
    orientation_file = sys.argv[2] if len(sys.argv) > 2 else default_orientation_file
    output_file = sys.argv[3] if len(sys.argv) > 3 else default_output_file

    print(f"\nProcessing centroid file: {centroid_file}")
    print(f"Processing orientation file: {orientation_file}")
    print(f"Output will be written to: {output_file}\n")

    # Initialization with updated file paths
    centroid = Centroid(centroid_file)
    tiltedorientation = TiltedOrientation(orientation_file)

    module_four_corners_database = {}

    for detid in centroid.data:
        x, y, z, moduleTypeInfo = centroid.getCentroid(detid)
        module = m.Module(detid, moduleTypeInfo)

        r = np.sqrt(x**2 + y**2)
        phi = np.arctan2(y, x)
        phi_ = np.arctan2(x, -y)
        is2S = module.moduleType()
        philen = 5
        etalen = 5 if is2S else 2.5

        if module.subdet() == 5: # Barrel
            drdz = get_drdz(tiltedorientation, detid, module.partnerDetId())
            four_corners = compute_four_corners_barrel(x, y, z, phi, phi_, etalen, philen, drdz, module.side())
        else: # Endcap
            four_corners = compute_four_corners_endcap(x, y, z, phi, phi_, etalen, philen)

        module_four_corners_database[detid] = four_corners

    with open(output_file, "w") as g:
        g.write(json.dumps(module_four_corners_database, indent=4))
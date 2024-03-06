import os
import sys
import json
import numpy as np
import pandas as pd

import Module as m
from compute_centroids import parse_module_type

# Constants for Barrel and Endcap subdet type
BARREL, ENDCAP = 5, 4
# Default value for slopes if undefined
DEFAULT_SLOPE = 123456789

def calculate_slope(dx, dy, dz):
    """
    Calculate the slopes drdz and dxdy based on the differences in x, y, and z corner coordinates.
    If slope is undefined, use the default slope value.
    
    Parameters:
    dx (float): Difference in x corner coordinates.
    dy (float): Difference in y corner coordinates.
    dz (float): Difference in z corner coordinates.
    
    Returns:
    tuple: A tuple containing the drdz_slope and dxdy_slope values.
    """
    dr = np.sqrt(dx**2 + dy**2)
    drdz_slope = dr / dz if dz != 0 else DEFAULT_SLOPE
    dxdy_slope = -dx / dy if dy != 0 else DEFAULT_SLOPE
    return drdz_slope, dxdy_slope

def process_corners(corners):
    """
    Use each sensor's corners to calculate and categorize drdz and dxdy slopes.
    
    Parameters:
    corners (dict): A dictionary with sensor DetID as keys and their corner coordinates as values.
    
    Returns:
    tuple: Two dictionaries containing the slopes for barrel and endcap sensors respectively.
    """
    barrel_slopes, endcap_slopes = {}, {}

    for detid, row in corners.items():
        # Grab first two corners of the sensor
        corner_1, corner_2 = row[0], row[1]

        # Compute dx, dy, dz for first two corner coordinates
        dx = corner_2[1] - corner_1[1]
        dy = corner_2[2] - corner_1[2]
        dz = corner_2[0] - corner_1[0]

        # Calculate drdz and dxdy slopes using the dx, dy, dz variables
        drdz_slope, dxdy_slope = calculate_slope(dx, dy, dz)

        # Initialize module object to get subdet, tilted and strip variables
        moduletype = parse_module_type(int(detid))
        module = m.Module(int(detid), int(moduletype))

        # Distinguishes barrel and endcap
        subdet = module.subdet()
        # Is the sensor tilted
        tilted = module.side() != 3
        # Is the sensor a strip sensor
        strip = module.moduleLayerType() == 1

        if not strip:
            continue

        # Save slopes to relevant dict if the sensor meets requirements
        slope_data = {'drdz_slope': drdz_slope, 'dxdy_slope': dxdy_slope}
        if subdet == BARREL and tilted:
            barrel_slopes[detid] = slope_data
        elif subdet == ENDCAP:
            endcap_slopes[detid] = slope_data

    return barrel_slopes, endcap_slopes

def save_slopes_to_file(slopes, output_path, centroids_df):
    """
    Saves the calculated slopes to a file.
    
    Parameters:
    slopes (dict): A dictionary with sensor DetID as keys and slopes as values.
    output_path (str): The file path to save the slope data.
    """
    with open(output_path, 'w') as file:
        for detid, values in slopes.items():
            # The endcap sensors will have drdz as the default always
            if values['drdz_slope'] != DEFAULT_SLOPE:
                line = f"{detid} {values['drdz_slope']} {values['dxdy_slope']}\n"
            else:
                centroid_phi = np.radians(centroids_df.loc[centroids_df['DetId/i'] == int(detid), ' phi_deg/D'].values[0])
                line = f"{detid} {values['dxdy_slope']} {centroid_phi}\n"
            file.write(line)

if __name__ == "__main__":
    # Default file paths
    default_sensor_corners_file = 'output/sensor_corners.txt'
    default_centroids_file = "data/DetId_sensors_list_OT800_IT615.csv"
    output_tilted_barrel_file = 'output/tilted_barrel_orientation.txt'
    output_endcap_file = 'output/endcap_orientation.txt'

    # Check for help flag
    if '-h' in sys.argv or '--help' in sys.argv:
        print("\nUsage: python compute_orientation.py [sensor_corners_file] [centroids_file] [output_tilted_barrel_file] [output_endcap_file]")
        print("\nOptions:")
        print(f"  sensor_corners_file        Path to the sensor corners file. Default is {default_sensor_corners_file}")
        print(f"  centroids_file             Path to the centroids file containing detId and phi values. Default is {default_centroids_file}")
        print(f"  output_tilted_barrel_file  Path to the output file for tilted barrel orientations. Default is {output_tilted_barrel_file}")
        print(f"  output_endcap_file         Path to the output file for endcap orientations. Default is {output_endcap_file}")
        sys.exit()

    # Determine file paths based on arguments provided
    sensor_corners_file = sys.argv[1] if len(sys.argv) > 1 else default_sensor_corners_file
    centroids_file = sys.argv[2] if len(sys.argv) > 2 else default_centroids_file
    output_tilted_barrel_file = sys.argv[3] if len(sys.argv) > 3 else output_tilted_barrel_file
    output_endcap_file = sys.argv[4] if len(sys.argv) > 4 else output_endcap_file

    # Ensure geometry file exists or create it
    if not os.path.exists(sensor_corners_file):
        from compute_geometry import main as make_geometry
        make_geometry()

    # Load sensor corners from file
    with open(sensor_corners_file, 'r') as file:
        corners = json.load(file)

    # Load centroids information for phi values
    centroids_df = pd.read_csv(default_centroids_file, usecols=['DetId/i', ' phi_deg/D'])

    # Process corners to calculate slopes
    barrel_slopes, endcap_slopes = process_corners(corners)

    # Save calculated slopes to files
    save_slopes_to_file(barrel_slopes, output_tilted_barrel_file, centroids_df)
    save_slopes_to_file(endcap_slopes, output_endcap_file, centroids_df)
import os
import sys
import json
import pandas as pd

# Import relevant functions from other scripts
from compute_corners import transform_sensor_corners, assign_corners_to_sensor
from compute_centroids import compute_centroids
from compute_orientation import process_corners, save_slopes_to_file

def compute_geometry(module_info_path, sensor_info_path, output_path_corners, output_path_centroid, output_path_tilted_barrel, output_path_endcap):
    """
    Computes the full geometry outputs including corners and centroids, and calculates orientations (slopes) for both barrel and endcap sensors.
    The results are saved to specified output files.
    
    Parameters:
    module_info_path (str): Path to the CSV file containing module information.
    sensor_info_path (str): Path to the CSV file containing sensor information.
    output_path_corners (str): File path for saving computed corners.
    output_path_centroid (str): File path for saving computed centroids.
    output_path_tilted_barrel (str): File path for saving calculated slopes for tilted barrel sensors.
    output_path_endcap (str): File path for saving calculated slopes for endcap sensors.
    
    This function integrates several steps including reading module and sensor information from CSV files,
    computing and assigning transformed corners to sensors, calculating centroids, and finally computing and saving orientations for barrel and endcap sensors.
    """
    # Read input module/sensor csv files
    module_csv = pd.read_csv(module_info_path)
    sensor_csv = pd.read_csv(sensor_info_path)

    # Compute corners and save to file
    module_csv['Transformed_Corners'] = module_csv.apply(transform_sensor_corners, axis=1)
    assigned_corners = assign_corners_to_sensor(module_csv, sensor_csv)
    with open(output_path_corners, 'w') as f:
        json.dump(assigned_corners, f, indent=4)

    # Compute centroids and save to file
    x, y, z, detid, moduleType = compute_centroids(sensor_info_path)
    with open(output_path_centroid, "w") as output:
        for i in range(len(x)):
            output.write(f"{detid[i]},{x[i]},{y[i]},{z[i]},{moduleType[i]}\n")

    # Compute barrel and endcap orientations (slopes) and save to files
    barrel_slopes, endcap_slopes = process_corners(assigned_corners)
    save_slopes_to_file(barrel_slopes, output_path_tilted_barrel, sensor_csv)
    save_slopes_to_file(endcap_slopes, output_path_endcap, sensor_csv)

    print(f"\nProcessed files: {module_info_path}, {sensor_info_path}")
    print(f"Output written to: {output_path_corners}, {output_path_centroid}, {output_path_tilted_barrel}, {output_path_endcap}\n")

def main():
    # Default file paths
    default_module_info_path = "data/module_info_OT800_IT615.csv"
    default_sensor_info_path = "data/DetId_sensors_list_OT800_IT615.csv"
    default_output_path_corners = "output/sensor_corners.txt"
    default_output_path_centroid = "output/sensor_centroids.txt"
    default_output_path_tilted_barrel = "output/tilted_barrel_orientation.txt"
    default_output_path_endcap = "output/endcap_orientation.txt"

    # Check for help flag
    if '-h' in sys.argv or '--help' in sys.argv:
        print("\nUsage: python compute_geometry.py [module_info_file] [sensor_info_file] [outputfile_corners] [outputfile_centroid] [outputfile_tilted_barrel] [outputfile_endcap]")
        print("\nOptions:")
        print(f"  module_info_file         Path to the module information CSV file. Default is {default_module_info_path}")
        print(f"  sensor_info_file         Path to the sensor information CSV file. Default is {default_sensor_info_path}")
        print(f"  outputfile_corners       Path for the corners output file. Default is {default_output_path_corners}")
        print(f"  outputfile_centroid      Path for the centroid output file. Default is {default_output_path_centroid}")
        print(f"  outputfile_tilted_barrel Path for the tilted barrel orientation output file. Default is {default_output_path_tilted_barrel}")
        print(f"  outputfile_endcap        Path for the endcap orientation output file. Default is {default_output_path_endcap}\n")
        sys.exit()

    # Determine input and output file paths based on arguments provided
    module_info_path = sys.argv[1] if len(sys.argv) > 1 else default_module_info_path
    sensor_info_path = sys.argv[2] if len(sys.argv) > 2 else default_sensor_info_path
    output_path_corners = sys.argv[3] if len(sys.argv) > 3 else default_output_path_corners
    output_path_centroid = sys.argv[4] if len(sys.argv) > 4 else default_output_path_centroid
    output_path_tilted_barrel = sys.argv[5] if len(sys.argv) > 5 else default_output_path_tilted_barrel
    output_path_endcap = sys.argv[6] if len(sys.argv) > 6 else default_output_path_endcap

    # Make output folder if it doesn't exist
    os.makedirs(os.path.dirname("output/"), exist_ok=True)

    # Compute geometry with specified file paths
    compute_geometry(module_info_path,
                     sensor_info_path,
                     output_path_corners,
                     output_path_centroid,
                     output_path_tilted_barrel,
                     output_path_endcap)

if __name__ == "__main__":
    main()
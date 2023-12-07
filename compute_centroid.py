import sys
import numpy as np
import pandas as pd

# Function to extract bits from detId
def extract_bits(value, start, end):
    mask = (1 << (end - start + 1)) - 1
    return (value >> start) & mask

def parse_module_type(det_id):
    """
    Get the numerical module type from the given detId.
    """
    # Check if the first digit of detId is '3' for inner tracker
    if str(det_id)[0] == '3':
        return -1

    # Define constants for subdetectors
    BARREL = 5
    ENDCAP = 4

    # Constants for module types
    PS = 23
    PSS = 24
    TwoS = 25

    # Parse subdet, layer, and ring
    subdet = extract_bits(det_id, 25, 27)
    layer = extract_bits(det_id, 18, 20) if subdet == ENDCAP else extract_bits(det_id, 20, 22)
    ring = extract_bits(det_id, 12, 15) if subdet == ENDCAP else 0

    # Determine module type based on subdet, layer, and ring
    if subdet == BARREL:
        return PS if layer <= 3 else TwoS
    elif subdet == ENDCAP:
        if layer <= 2:
            return PS if ring <= 10 else TwoS
        else:
            return PS if ring <= 7 else TwoS
    else:
        raise ValueError("Invalid subdet value")

def process_csv(file_path):
    df_parsed = pd.read_csv(file_path)

    # Extracting rho, phi, Z, and detid from the dataframe
    rho = df_parsed[' sensorCenterRho_mm/D']
    phi_deg = df_parsed[' phi_deg/D']
    z = df_parsed[' sensorCenterZ_mm/D']
    detid = df_parsed['DetId/i']

    # Extract moduleType for each detid
    moduleType = np.array([parse_module_type(d) for d in detid])

    # Remove sensors from inner tracker
    itmask = moduleType != -1
    rho, phi_deg, z, detid, moduleType = rho[itmask], phi_deg[itmask], z[itmask], detid[itmask], moduleType[itmask]

    # Converting phi from degrees to radians and calculating X and Y from rho, phi
    phi_rad = np.radians(phi_deg)
    x = rho * np.cos(phi_rad)
    y = rho * np.sin(phi_rad)

    # Account for scaling difference used in excel sheet.
    x, y, z = x/10, y/10, z/10

    return x, y, z, detid, moduleType

if __name__ == "__main__":
    # Default file paths
    default_input_path = "data/DetId_sensors_list.csv"
    default_output_path = "data/centroid.txt"

    # Check for help flag
    if '-h' in sys.argv or '--help' in sys.argv:
        print("\nUsage: python compute_centroid.py [inputfile] [outputfile]")
        print("\nOptions:")
        print(f"  inputfile   Path to the input CSV file. Default is {default_input_path}")
        print(f"  outputfile  Path for the output file. Default is {default_output_path}\n")
        sys.exit()

    # Determine input and output file paths based on arguments provided
    input_path = sys.argv[1] if len(sys.argv) > 1 else default_input_path
    output_path = sys.argv[2] if len(sys.argv) > 2 else default_output_path

    # Process CSV file
    x, y, z, detid, moduleType = process_csv(input_path)

    # Write the data to the specified output file
    with open(output_path, "w") as output:
        for i in range(len(x)):
            output.write(f"{detid[i]},{x[i]},{y[i]},{z[i]},{moduleType[i]}\n")

    print(f"\nProcessed file: {input_path}")
    print(f"Output written to: {output_path}\n")
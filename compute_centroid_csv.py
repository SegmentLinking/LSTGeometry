import sys
import numpy as np
import pandas as pd

def extract_bits(value, start, end):
    """
    Extracts a specific range of bits from an integer value.

    Parameters:
    value (int): The integer from which bits will be extracted.
    start (int): The starting bit position (inclusive).
    end (int): The ending bit position (inclusive).

    Returns:
    int: The extracted bits as an integer.
    """
    mask = (1 << (end - start + 1)) - 1
    return (value >> start) & mask

def parse_module_type(det_id):
    """
    Determines the module type of a sensor based on its detector ID.

    Parameters:
    det_id (int): The detector ID of the sensor.

    Returns:
    int: The numerical module type. 
         -1 for inner tracker modules,
         23 (PSP), 24 (PSS), or 25 (TwoS) for different module types.
    
    Note: 
    For more information, see
    https://github.com/cms-sw/cmssw/tree/master/Geometry/TrackerGeometryBuilder
    """
    # Check if the first digit of detId is '3' for inner tracker
    if str(det_id)[0] == '3':
        return -1

    # Define constants for subdetectors
    BARREL = 5
    ENDCAP = 4

    # Constants for module types
    PSP = 23
    PSS = 24
    TwoS = 25

    # Parse subdet, layer, and ring
    subdet = extract_bits(det_id, 25, 27)
    layer = extract_bits(det_id, 20, 22) if subdet == BARREL else extract_bits(det_id, 18, 20)
    ring = extract_bits(det_id, 12, 15) if subdet == ENDCAP else 0

    # Determine module type based on subdet, layer, and ring
    is_even_det_id = det_id % 2 == 0
    if subdet == BARREL:
        if layer <= 3:
            return PSS if is_even_det_id else PSP
        else:
            return TwoS
    elif subdet == ENDCAP:
        if layer <= 2:
            return PSS if ring <= 10 and is_even_det_id else (PSP if ring <= 10 else TwoS)
        else:
            return PSS if ring <= 7 and is_even_det_id else (PSP if ring <= 7 else TwoS)
    else:
        raise ValueError("Invalid subdet value")

def process_csv(file_path):
    """
    Use CSV file containing sensor data (rho, phi, Z, detid) to return sensor centroid coordinates (X, Y, Z).

    Parameters:
    file_path (str): The path to the sensor CSV file to be processed.

    Returns:
    tuple: A tuple containing the processed data arrays (x, y, z, detid, moduleType).

    Note:
    Excel spreadsheet values are given in mm, but later scripts expect this scaled by 1/10.
    """
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
        print("\nUsage: python compute_centroid_csv.py [inputfile] [outputfile]")
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
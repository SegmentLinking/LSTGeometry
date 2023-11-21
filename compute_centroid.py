import sys
import numpy as np

def parse_hits(file_path):
    """
    Parses the input txt file to extract hits data.

    Parameters:
    file_path (str): Path to the input file containing hits data.

    Returns:
    dict: A dictionary containing NumPy arrays for x, y, z coordinates, detId, and moduleType.
    """
    x, y, z, detIds, moduleTypes = [], [], [], [], []
    uniqueDetIds = set()

    with open(file_path) as f:
        lines = f.readlines()

    for line in lines:
        if "detId" not in line:
            continue
        ls = line.split()

        # Append values for this line to relevant arrays
        x.append(float(ls[1]))  # x coordinate
        y.append(float(ls[3]))  # y coordinate
        z.append(float(ls[5]))  # z coordinate
        detIds.append(int(ls[7]))
        moduleTypes.append(int(ls[9]))

    return {'x': np.array(x, dtype=np.float32), 
            'y': np.array(y, dtype=np.float32), 
            'z': np.array(z, dtype=np.float32), 
            'detId': np.array(detIds, dtype=np.uint32), 
            'moduleType': np.array(moduleTypes, dtype=np.uint32)}

def calculate_centroids(hit_data):
    """
    Calculates centroids for each unique detector ID (detId) using the provided hit data.

    Parameters:
    hit_data (dict): A dictionary containing NumPy arrays for x, y, z coordinates, and detId.

    Returns:
    dict: A dictionary with centroids (mean positions) for each unique detId.
    """
    unique_detIds = np.unique(hit_data['detId'])

    # Initialize centroids
    centroids = {}

    # Calculate weighted sums for x, y, z coordinates
    sums_x = np.bincount(hit_data['detId'], weights=hit_data['x'])
    sums_y = np.bincount(hit_data['detId'], weights=hit_data['y'])
    sums_z = np.bincount(hit_data['detId'], weights=hit_data['z'])

    # Count occurrences of each detId
    counts = np.bincount(hit_data['detId'])

    # Calculate centroids
    for detId in unique_detIds:
        centroids[detId] = np.array([
            sums_x[detId] / counts[detId], 
            sums_y[detId] / counts[detId], 
            sums_z[detId] / counts[detId]
        ])

    return centroids

def calculate_module_types(hit_data):
    """
    Extracts module types for each unique detector ID (detId) from the provided hit data.

    Parameters:
    hit_data (dict): A dictionary containing NumPy arrays for detId and moduleType.

    Returns:
    dict: A dictionary mapping each unique detId to its module type.
    """
    # Extract unique detector IDs and their first occurrence indices
    unique_detIds, first_indices = np.unique(hit_data['detId'], return_index=True)

    # Map each unique detId to its corresponding moduleType using the first occurrence
    moduleTypes = {hit_data['detId'][index]: hit_data['moduleType'][index] for index in first_indices}

    return moduleTypes

if __name__ == "__main__":
    # Default file paths
    default_input_path = "./hits.txt"
    default_output_path = "data/centroid.txt"

    # Check for help flag
    if '-h' in sys.argv or '--help' in sys.argv:
        print("\nUsage: python script_name.py [inputfile] [outputfile]")
        print("\nOptions:")
        print(f"  inputfile   Path to the input TXT file containing hits data. Default is {default_input_path}")
        print(f"  outputfile  Path for the output file. Default is {default_output_path}\n")
        sys.exit()

    # Determine input and output file paths based on arguments provided
    input_path = sys.argv[1] if len(sys.argv) > 1 else default_input_path
    output_path = sys.argv[2] if len(sys.argv) > 2 else default_output_path

    print(f"\nProcessing file: {input_path}")
    print(f"Output will be written to: {output_path}\n")

    # Get hit info from input txt file
    hit_data = parse_hits(input_path)
    
    # Calculate centroids and moduleType dictionary
    centroids = calculate_centroids(hit_data)
    moduleType = calculate_module_types(hit_data)

    # Write the centroid data and moduletype to the specified output file
    with open(output_path, "w") as output:
        for detid, centroid in centroids.items():
            output.write(f"{detid},{centroid[0]},{centroid[1]},{centroid[2]},{moduleType[detid]}\n")
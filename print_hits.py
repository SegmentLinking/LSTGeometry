import sys

def process_hits(input_filename, output_filename=None):
    """
    Process hit data from a ROOT file and either write to an output file or collect to return.

    Parameters:
    - input_filename: str, the path to the ROOT file.
    - output_filename: Optional[str], the path to the output text file.

    Returns:
    - A list of dictionaries with hit data if output_filename is None, otherwise None.
    """
    import uproot
    import numpy as np

    # Prepare branches to read
    branches = ["ph2_x", "ph2_y", "ph2_z", "ph2_detId", "ph2_moduleType"]

    # Define a function to process a single batch of hits
    def process_batch(batch, output_file=None):
        hit_data_batch = []
        for i in range(len(batch["ph2_x"])):
            for hit in range(len(batch["ph2_x"][i])):
                hit_info = {
                    "x": batch["ph2_x"][i][hit],
                    "y": batch["ph2_y"][i][hit],
                    "z": batch["ph2_z"][i][hit],
                    "detId": batch["ph2_detId"][i][hit],
                    "moduleType": batch["ph2_moduleType"][i][hit]
                }
                if output_file:
                    output_file.write(f"x: {hit_info['x']} y: {hit_info['y']} z: {hit_info['z']} " +
                                      f"detId: {hit_info['detId']} moduleType: {hit_info['moduleType']}\n")
                else:
                    hit_data_batch.append(hit_info)
        return hit_data_batch

    # Process the file and write to output file or collect to return
    with uproot.open(input_filename) as file:
        tree = file["tree;11"]
        if output_filename:
            with open(output_filename, "w") as output_file:
                for batch in tree.iterate(branches, library="np"):
                    process_batch(batch, output_file)
        else:
            hit_data = []
            for batch in tree.iterate(branches, library="np"):
                hit_data.extend(process_batch(batch))
            return hit_data

if __name__ == "__main__":
    # Check command-line arguments
    if len(sys.argv) < 2:
        print("Usage: python script.py <inputfile> [outputfile]")
        sys.exit(1)

    # Input filename is required
    input_file = sys.argv[1]

    # Output filename is optional
    output_file = sys.argv[2] if len(sys.argv) > 2 else "hits.txt"

    process_hits(input_file, output_file)

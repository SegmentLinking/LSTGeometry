import sys
import uproot
import numpy as np

def process_hits(input_filename, output_filename=None, nevents=100):
    """
    Process hit data from a ROOT file and either write to an output file or collect to return.

    Parameters:
    - input_filename: str, the path to the ROOT file.
    - output_filename: Optional[str], the path to the output text file.
    - nevents: int, the maximum number of events to process from the ROOT file.

    Returns:
    - A dictionary of arrays with hit data if output_filename is None, otherwise None.
    """
    # Prepare branches to read
    branches = ["ph2_x", "ph2_y", "ph2_z", "ph2_detId", "ph2_moduleType"]

    # Define a function to process a single batch of hits
    def process_batch(batch, output_file=None):
        if output_file:
            for x, y, z, detId, moduleType in zip(batch["ph2_x"], batch["ph2_y"], batch["ph2_z"], batch["ph2_detId"], batch["ph2_moduleType"]):
                for hit in range(len(x)):
                    output_file.write(f"x: {x[hit]} y: {y[hit]} z: {z[hit]} detId: {detId[hit]} moduleType: {moduleType[hit]}\n")
        else:
            batch_data = {key: [] for key in branches}
            for key in branches:
                batch_data[key] = np.concatenate(batch[key])
            return batch_data

    # Process the file and write to output file or collect to return
    with uproot.open(input_filename) as file:
        try:
            tree = file["trackingNtuple/tree"]
        except KeyError:
            # If trackingNtuple/tree is not found, try to find another tree
            available_trees = list(file.keys())
            if not available_trees:
                raise RuntimeError("No trees found in the ROOT file.")
            tree = file[available_trees[0]]  # Try the first available tree

        if output_filename:
            with open(output_filename, "w") as output_file:
                for batch in tree.iterate(branches, library="np", entry_stop=nevents):
                    process_batch(batch, output_file)
        else:
            hit_data = {key: [] for key in branches}
            for batch in tree.iterate(branches, library="np", entry_stop=nevents):
                batch_data = process_batch(batch)
                for key in hit_data:
                    hit_data[key].extend(batch_data[key])
            return {k: np.array(v) for k, v in hit_data.items()}

if __name__ == "__main__":
    default_input_file = "/data2/segmentlinking/CMSSW_12_5_0_pre3/RelValTTbar_14TeV_CMSSW_12_5_0_pre3/event_1000.root"

    # Check for help flag
    if '-h' in sys.argv or '--help' in sys.argv:
        print("\nUsage: python print_hits.py [inputfile] [outputfile] [nevents]")
        print("\nOptions:")
        print(f"  inputfile   The path to the ROOT file to process. Default is {default_input_file}")
        print("  outputfile  The path to the output text file. Default is hits.txt")
        print("  nevents     The maximum number of events to process. Default is 100\n")
        sys.exit()

    # Determine input/output filenames and number of events based on arguments provided
    input_file = sys.argv[1] if len(sys.argv) > 1 else default_input_file
    output_file = sys.argv[2] if len(sys.argv) > 2 else "hits.txt"
    max_events = int(sys.argv[3]) if len(sys.argv) > 3 else 100

    print(f"\nProcessing file: {input_file}")
    print(f"Output will be written to: {output_file}")
    print(f"Maximum number of events to process: {max_events}\n")

    process_hits(input_file, output_file, max_events)

import sys
import tqdm
import numpy as np
import Module as m
from scipy.optimize import curve_fit

def parse_and_store_hits(filename):
    """
    Parses the hits from a file and stores them in a dictionary.

    Parameters:
    - filename: The path to the input text file.

    Returns:
    - hits: A dictionary containing hits grouped by detector ID.
    - moduleType: A dictionary mapping detector IDs to module types.
    """
    hits = {}
    moduleType = {}

    with open(filename) as f:
        for line in f:
            if "detId" not in line:
                continue
            ls = line.split()
            detid = ls[7]
            moduleType.setdefault(detid, ls[9])
            hit = (float(ls[1]), float(ls[3]), float(ls[5]))
            hits.setdefault(detid, []).append(hit)

    return hits, moduleType

def compute_hit_groups(detid, hits):
    """
    Groups hits based on their z-coordinates into low-z and high-z groups.

    Parameters:
    - detid: The detector ID.
    - hits: List of hits for the given detector ID.

    Returns:
    - low_z_hits: List of hits with lower z-coordinates.
    - high_z_hits: List of hits with higher z-coordinates.
    - x_low_z: List of x-coordinates for low-z hits.
    - x_high_z: List of x-coordinates for high-z hits.
    - sorted_abs_z: Sorted list of unique absolute z-coordinates.
    """
    sorted_abs_z = sorted(set(abs(hit[2]) for hit in hits[detid]))
    low_z_hits, high_z_hits = [], []
    x_low_z, x_high_z = [], []
    for hit in hits[detid]:
        az = abs(hit[2])
        if az == sorted_abs_z[0]:
            low_z_hits.append(hit)
            x_low_z.append(hit[0])
        else:
            high_z_hits.append(hit)
            x_high_z.append(hit[0])
    return low_z_hits, high_z_hits, x_low_z, x_high_z, sorted_abs_z

def linear_fit_func(x, a, b):
    return a * x + b

def fit_hits(hits):
    """
    Fits a polynomial to the provided hits and returns the fit parameters.
    Handles special case where all x-coordinates are the same (vertical line).

    Parameters:
    - hits: List of hits to fit.

    Returns:
    - y_intercept: The y-intercept from the fit.
    - slope: The slope from the fit, or a special value for vertical lines.
    """
    if len(hits) < 2:
        return -999, -999

    x = np.array([hit[0] for hit in hits])
    y = np.array([hit[1] for hit in hits])

    # Check if all x-coordinates are the same (vertical line)
    if np.all(x == x[0]):
        return -999, -999

    try:
        params, _ = curve_fit(linear_fit_func, x, y)
        return params[1], params[0]  # y_intercept, slope
    except RuntimeError:
        return -999, -999

if __name__ == "__main__":
    # Default file paths
    default_input_file = "hits.txt"
    default_output_file = "data/tilted_orientation.txt"

    # Check for help flag
    if '-h' in sys.argv or '--help' in sys.argv:
        print("\nUsage: python fit_tilted_module_orientations.py [inputfile] [outputfile]")
        print("\nOptions:")
        print(f"  inputfile   Path to the input TXT file containing hits data. Default is {default_input_file}")
        print(f"  outputfile  Path for the output file. Default is {default_output_file}\n")
        sys.exit()

    # Determine input and output file paths based on arguments provided
    input_filename = sys.argv[1] if len(sys.argv) > 1 else default_input_file
    output_filename = sys.argv[2] if len(sys.argv) > 2 else default_output_file

    print(f"\nProcessing file: {input_filename}")
    print(f"Output will be written to: {output_filename}\n")

    # Pull hit data, moduletype from input file
    hits, moduleType = parse_and_store_hits(input_filename)

    with open(output_filename, "w") as output:
        output.write("# detid drdz xy-slope\n")

        # Loop through detid's and calculate orientation of each.
        for detid in tqdm.tqdm(hits):
            mod = m.Module(int(detid), int(moduleType[detid]))

            istilt = (mod.side() == 1 or mod.side() == 2) and mod.subdet() == 5
            if not istilt:
                continue

            isstrip = mod.moduleLayerType() == 1
            if not isstrip:
                continue

            low_z_hits, high_z_hits, x_low_z, x_high_z, sorted_abs_z = compute_hit_groups(detid, hits)

            # Get y-intercepts and slopes from fit
            yl, sl = fit_hits(low_z_hits)
            yh, sh = fit_hits(high_z_hits)

            if yl != -999 and yh != -999:
                dydz = abs(yh - yl) / np.sqrt(sl**2 + 1) / abs(sorted_abs_z[0] - sorted_abs_z[1])
                output.write(f"{detid} {dydz} {sl}\n")
            else:
                dxdz = abs(list(set(x_low_z))[0] - list(set(x_high_z))[0]) / abs(sorted_abs_z[0] - sorted_abs_z[1])
                output.write(f"{detid} {dxdz} 123456789\n")
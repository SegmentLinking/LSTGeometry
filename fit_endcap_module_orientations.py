import sys
import tqdm
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
    Computes two groups of hits based on the average value of radial distance squared.

    Parameters:
    - detid: The detector ID.
    - hits: The list of hits for the given detector ID.

    Returns:
    - low_hits: Hits with lower radial distance squared.
    - high_hits: Hits with higher radial distance squared.
    """
    n = len(hits[detid])
    sumr2s = sum(hit[0]**2 + hit[1]**2 for hit in hits[detid])
    avgr2s = sumr2s / n

    low_hits = [hit for hit in hits[detid] if (hit[0]**2 + hit[1]**2) < avgr2s]
    high_hits = [hit for hit in hits[detid] if (hit[0]**2 + hit[1]**2) >= avgr2s]

    return low_hits, high_hits

def linear_fit_func(x, a, b):
    return a * x + b

def fit_hits(hits):
    """
    Fits a polynomial to the provided hits and returns the fit parameters.

    Parameters:
    - hits: List of hits to fit.

    Returns:
    - y_intercept: The y-intercept from the fit.
    - slope: The slope from the fit.
    """
    if len(hits) < 2:
        return -999, -999

    # Extract x and y coordinates
    x = [hit[0] for hit in hits]
    y = [hit[1] for hit in hits]

    try:
        # Fit the data
        params, _ = curve_fit(linear_fit_func, x, y)
        return params[1], params[0]  # y_intercept, slope
    except RuntimeError:
        # Fit failed
        return -999, -999

if __name__ == "__main__":
    # Default file paths
    default_input_file = "hits.txt"
    default_output_file = "data/endcap_orientation.txt"

    # Check for help flag
    if '-h' in sys.argv or '--help' in sys.argv:
        print("\nUsage: python fit_endcap_module_orientations.py [inputfile] [outputfile]")
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
        output.write("# detid average_r2s y_intercept_low_hits slope_low_hits y_intercept_high_hits slope_high_hits\n")

        # Loop through detid's and calculate orientation of each.
        for detid in tqdm.tqdm(hits):
            mod = m.Module(int(detid), int(moduleType[detid]))

            isendcap = mod.subdet() == 4
            if not isendcap:
                continue

            isstrip = mod.moduleLayerType() == 1
            if not isstrip:
                continue

            low_hits, high_hits = compute_hit_groups(detid, hits)

            # Get y-intercepts and slopes from fit
            yl, sl = fit_hits(low_hits)
            yh, sh = fit_hits(high_hits)

            # Adjustments for missing fits
            if sl == -999 and sh != -999:
                sl, yl = sh, yh
            elif sl != -999 and sh == -999:
                sh, yh = sl, yl

            output.write(f"{detid} {sum(hit[0]**2 + hit[1]**2 for hit in hits[detid]) / len(hits[detid])} {yl} {sl} {yh} {sh}\n")
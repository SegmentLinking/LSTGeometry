import ROOT as r
import json
import os
import sys
from tqdm import tqdm
from Module import Module
from Centroid import Centroid
import math
from DetectorGeometry import DetectorGeometry

# See Constants.py for definitions.
from Constants import PTTHRESH

def printPixelMap(centroid, det_geom, output_dir="output/pixelmap/"):
    """
    To print out pixel maps
    """

    # Define the binning of "super bins"
    neta = 25.
    nphi = 72.
    nz = 25.
    # pt_bounds = [0.9, 2.0, 4.0, 10., 50.]
    pt_bounds = [PTTHRESH, 2.0, 10000.]

    # Grand map object that will hold the mapping of the pixel map
    # maps[(ipt, ieta, iphi, iz)][(layer, subdet)] = [] # for both positive and negative
    # maps[(ipt, ieta, iphi, iz, 1)][(layer, subdet)] = [] # for positive oly
    # maps[(ipt, ieta, iphi, iz,-1)][(layer, subdet)] = [] # for negative oly
    maps = {}

    # Initialize empty lists for the pixel map 
    for ipt in range(len(pt_bounds)-1):
        for ieta in range(int(neta)):
            for iphi in range(int(nphi)):
                for iz in range(int(nz)):

                    # Maps without split by charge
                    maps[(ipt, ieta, iphi, iz)] = {}

                    # Maps with split by charge (positive)
                    maps[(ipt, ieta, iphi, iz, 1)] = {}

                    # Maps with split by charge (negative)
                    maps[(ipt, ieta, iphi, iz, -1)] = {}

                    # The maps will then be split by (layer, subdet)
                    for layer in [1, 2]: #[1, 2, 3, 4, 5, 6]:
                        for subdet in [4, 5]:

                            # Maps without split by charge
                            maps[(ipt, ieta, iphi, iz)][(layer, subdet)] = []

                            # Maps split by charge
                            maps[(ipt, ieta, iphi, iz, 1)][(layer, subdet)] = []

                            # Maps split by charge
                            maps[(ipt, ieta, iphi, iz, -1)][(layer, subdet)] = []

    # Loop over the detids and for each detid compute which superbins it is connected to
    for detid in tqdm(det_geom.getDetIds()):
        x, y, z, moduleTypeInfo = centroid.getCentroid(detid)

        # Parse the layer and subdet
        module = Module(detid, moduleTypeInfo)
        layer = module.layer()
        if layer > 2: continue
        subdet = module.subdet()

        # Skip if the module is not PS module and is not lower module
        if module.isLower() != 1 or module.moduleType() != 0:
            continue

        # For this module, now compute which super bins they belong to
        # To compute which super bins it belongs to, one needs to provide at least pt and z window to compute compatible eta and phi range
        # So we have a loop in pt and Z
        for ipt in range(len(pt_bounds)-1):
            for iz in range(int(nz)):

                # The zmin, zmax of consideration
                zmin = -30 + iz * (60. / nz)
                zmax = -30 + (iz + 1) * (60. / nz)

                zmin -= 0.05
                zmin += 0.05

                # The ptmin, ptmax of consideration
                pt_lo = pt_bounds[ipt]
                pt_hi = pt_bounds[ipt+1]

                # Based on the zmin and zmax, and the detid we can get eta min and eta max of compatible range
                etamin = det_geom.getCompatibleEtaRange(detid, zmin, zmax)[0]
                etamax = det_geom.getCompatibleEtaRange(detid, zmin, zmax)[1]

                etamin -= 0.05
                etamax += 0.05

                if layer==2 and subdet==4:
                    if etamax < 2.3: continue
                    if etamin < 2.3: etamin = 2.3

                # Compute the indices of the compatible eta range
                ietamin = int((etamin + 2.6) / (5.2 / neta))
                ietamax = int((etamax + 2.6) / (5.2 / neta))

                # Since the etas are restricted to 2.6, need to chop it off if it's out of bounds
                # prelim_etabins = range(ietamin-1, ietamax+2) # add -1 and +1 to cover some inefficiencies
                prelim_etabins = range(ietamin, ietamax+1) # add -1 and +1 to cover some inefficiencies
                etabins = []
                for ieta in prelim_etabins:
                    if ieta >= 0 and ieta < neta:
                        etabins.append(ieta)

                # Now compute the ranges of iphi min and max for given pt ranges
                iphimin_pos = int((det_geom.getCompatiblePhiRange(detid, pt_lo, pt_hi)[0][0] + math.pi) / (2*math.pi / nphi))
                iphimax_pos = int((det_geom.getCompatiblePhiRange(detid, pt_lo, pt_hi)[0][1] + math.pi) / (2*math.pi / nphi))
                iphimin_neg = int((det_geom.getCompatiblePhiRange(detid, pt_lo, pt_hi)[1][0] + math.pi) / (2*math.pi / nphi))
                iphimax_neg = int((det_geom.getCompatiblePhiRange(detid, pt_lo, pt_hi)[1][1] + math.pi) / (2*math.pi / nphi))

                # if the range is crossing the -pi v. pi boundary special care is needed
                if iphimin_pos <= iphimax_pos:
                    phibins_pos = list(range(iphimin_pos, iphimax_pos))
                else:
                    phibins_pos = list(range(0, iphimax_pos)) + list(range(iphimin_pos, int(nphi)))
                if iphimin_neg <= iphimax_neg:
                    phibins_neg = list(range(iphimin_neg, iphimax_neg))
                else:
                    phibins_neg = list(range(0, iphimax_neg)) + list(range(iphimin_neg, int(nphi)))

                # Now we have a list of (ipt, ieta, iphi, iz)'s that are compatible so we fill the map
                for ieta in etabins:
                    for iphi in phibins_pos:
                        maps[(ipt, ieta, iphi, iz)][(layer, subdet)].append(detid)
                    for iphi in phibins_neg:
                        maps[(ipt, ieta, iphi, iz)][(layer, subdet)].append(detid)
                    for iphi in phibins_pos:
                        maps[(ipt, ieta, iphi, iz, 1)][(layer, subdet)].append(detid)
                    for iphi in phibins_neg:
                        maps[(ipt, ieta, iphi, iz, -1)][(layer, subdet)].append(detid)

    # Grand map txt file that will hold everything regardless of the layers
    g = open(output_dir + "pLS_map_ElCheapo.txt", "w")
    g_pos = open(output_dir + "pLS_map_pos_ElCheapo.txt", "w")
    g_neg = open(output_dir + "pLS_map_neg_ElCheapo.txt", "w")

    # pixel maps split by layers
    fs = {}
    for layer in [1, 2]: #[1, 2, 3, 4, 5, 6]:
        for subdet in [4, 5]:
            fs[(layer, subdet)] = open(output_dir + "pLS_map_layer{}_subdet{}.txt".format(layer, subdet), "w")
            fs[(layer, subdet, 1)] = open(output_dir + "pLS_map_pos_layer{}_subdet{}.txt".format(layer, subdet), "w")
            fs[(layer, subdet, -1)] = open(output_dir + "pLS_map_neg_layer{}_subdet{}.txt".format(layer, subdet), "w")

    # Loop over the super bins
    for ipt in range(len(pt_bounds)-1):
        for ieta in range(int(neta)):
            for iphi in range(int(nphi)):
                for iz in range(int(nz)):

                    # Compute the superbin index the phase-space belongs to
                    isuperbin = (ipt * nphi * neta * nz) + (ieta * nphi * nz) + nz * iphi + iz

                    # Temporary list to aggregate the detids
                    all_detids = []
                    all_pos_detids = []
                    all_neg_detids = []

                    # Loop over the possible layer and subdets
                    for layer in [1, 2]: #[1, 2, 3, 4, 5, 6]:
                        for subdet in [4, 5]:

                            # Obtain unique set of detids
                            maps[(ipt, ieta, iphi, iz)][(layer, subdet)] = list(set(maps[(ipt, ieta, iphi, iz)][(layer, subdet)]))
                            maps[(ipt, ieta, iphi, iz, 1)][(layer, subdet)] = list(set(maps[(ipt, ieta, iphi, iz, 1)][(layer, subdet)]))
                            maps[(ipt, ieta, iphi, iz, -1)][(layer, subdet)] = list(set(maps[(ipt, ieta, iphi, iz, -1)][(layer, subdet)]))

                            # Write out to the individual map files
                            fs[(layer, subdet)].write("{} {} {}\n".format(int(isuperbin), len(maps[(ipt, ieta, iphi, iz)][(layer, subdet)]), " ".join([ str(x) for x in maps[(ipt, ieta, iphi, iz)][(layer, subdet)] ])))
                            fs[(layer, subdet, 1)].write("{} {} {}\n".format(int(isuperbin), len(maps[(ipt, ieta, iphi, iz, 1)][(layer, subdet)]), " ".join([ str(x) for x in maps[(ipt, ieta, iphi, iz, 1)][(layer, subdet)] ])))
                            fs[(layer, subdet, -1)].write("{} {} {}\n".format(int(isuperbin), len(maps[(ipt, ieta, iphi, iz, -1)][(layer, subdet)]), " ".join([ str(x) for x in maps[(ipt, ieta, iphi, iz, -1)][(layer, subdet)] ])))

                            # Aggregate all the detids for a given superbin in order to later write out to the grand map txt file that holds all detids regardless of layers
                            all_detids += maps[(ipt, ieta, iphi, iz)][(layer, subdet)]
                            all_pos_detids += maps[(ipt, ieta, iphi, iz, 1)][(layer, subdet)]
                            all_neg_detids += maps[(ipt, ieta, iphi, iz, -1)][(layer, subdet)]

                    maps[(ipt, ieta, iphi, iz)]["all"] = list(set(all_detids))
                    maps[(ipt, ieta, iphi, iz, 1)]["all"] = list(set(all_pos_detids))
                    maps[(ipt, ieta, iphi, iz,-1)]["all"] = list(set(all_neg_detids))

                    # Now write the entire detid for the entire outer tracker regardless of layers for the superbin at hand
                    g.write("{} {} {}\n".format(int(isuperbin), len(maps[(ipt, ieta, iphi, iz)]["all"]), " ".join([ str(x) for x in maps[(ipt, ieta, iphi, iz)]["all"] ])))
                    g_pos.write("{} {} {}\n".format(int(isuperbin), len(maps[(ipt, ieta, iphi, iz, 1)]["all"]), " ".join([ str(x) for x in maps[(ipt, ieta, iphi, iz, 1)]["all"] ])))
                    g_neg.write("{} {} {}\n".format(int(isuperbin), len(maps[(ipt, ieta, iphi, iz,-1)]["all"]), " ".join([ str(x) for x in maps[(ipt, ieta, iphi, iz,-1)]["all"] ])))

    # Declaring histograms for mapping multiplicities
    ofile = r.TFile(output_dir + "pLS_map.root", "recreate")
    nconns = {}
    for ipt in range(len(pt_bounds)-1):
        for iz in range(int(nz)):
            nconns[(ipt, iz)] = r.TH2F("pLS_map_ElCheapo_ipt{}_iz{}".format(ipt, iz), "", int(nphi), -math.pi, math.pi, int(neta), -2.6, 2.6)
            nconns[(ipt, iz, 1)] = r.TH2F("pLS_map_pos_ElCheapo_ipt{}_iz{}".format(ipt, iz), "", int(nphi), -math.pi, math.pi, int(neta), -2.6, 2.6)
            nconns[(ipt, iz,-1)] = r.TH2F("pLS_map_neg_ElCheapo_ipt{}_iz{}".format(ipt, iz), "", int(nphi), -math.pi, math.pi, int(neta), -2.6, 2.6)
            for layer in [1, 2]: #[1, 2, 3, 4, 5, 6]:
                for subdet in [4, 5]:
                    nconns[(ipt, iz, layer, subdet)] = r.TH2F("pLS_map_layer{}_subdet{}_ipt{}_iz{}".format(layer, subdet, ipt, iz), "", int(nphi), -math.pi, math.pi, int(neta), -2.6, 2.6)
                    nconns[(ipt, iz, layer, subdet, 1)] = r.TH2F("pLS_map_pos_layer{}_subdet{}_ipt{}_iz{}".format(layer, subdet, ipt, iz), "", int(nphi), -math.pi, math.pi, int(neta), -2.6, 2.6)
                    nconns[(ipt, iz, layer, subdet,-1)] = r.TH2F("pLS_map_neg_layer{}_subdet{}_ipt{}_iz{}".format(layer, subdet, ipt, iz), "", int(nphi), -math.pi, math.pi, int(neta), -2.6, 2.6)

    # Now filling the histograms
    for ipt in range(len(pt_bounds)-1):
        for iz in range(int(nz)):
            for ieta in range(int(neta)):
                for iphi in range(int(nphi)):
                    nconns[(ipt, iz)].SetBinContent(iphi + 1, ieta + 1, len(maps[(ipt, ieta, iphi, iz)]["all"]))
                    nconns[(ipt, iz, 1)].SetBinContent(iphi + 1, ieta + 1, len(maps[(ipt, ieta, iphi, iz, 1)]["all"]))
                    nconns[(ipt, iz,-1)].SetBinContent(iphi + 1, ieta + 1, len(maps[(ipt, ieta, iphi, iz,-1)]["all"]))
            for layer in [1, 2]: #[1, 2, 3, 4, 5, 6]:
                for subdet in [4, 5]:
                    for ieta in range(int(neta)):
                        for iphi in range(int(nphi)):
                            nconns[(ipt, iz, layer, subdet)].SetBinContent(iphi + 1, ieta + 1, len(maps[(ipt, ieta, iphi, iz)][(layer, subdet)]))
                            nconns[(ipt, iz, layer, subdet, 1)].SetBinContent(iphi + 1, ieta + 1, len(maps[(ipt, ieta, iphi, iz, 1)][(layer, subdet)]))
                            nconns[(ipt, iz, layer, subdet,-1)].SetBinContent(iphi + 1, ieta + 1, len(maps[(ipt, ieta, iphi, iz,-1)][(layer, subdet)]))

    # Write the Histograms
    ofile.cd()
    for ipt in range(len(pt_bounds)-1):
        for iz in range(int(nz)):
            nconns[(ipt, iz)].Write()
            nconns[(ipt, iz, 1)].Write()
            nconns[(ipt, iz,-1)].Write()
            for layer in [1, 2]: #[1, 2, 3, 4, 5, 6]:
                for subdet in [4, 5]:
                    nconns[(ipt, iz, layer, subdet)].Write()
                    nconns[(ipt, iz, layer, subdet, 1)].Write()
                    nconns[(ipt, iz, layer, subdet,-1)].Write()

if __name__ == "__main__":
    # Default file paths
    default_centroid_file = "output/sensor_centroids.txt"
    default_geom_file = "output/sensor_corners.txt"
    default_average_radius_file = "data/average_r_OT800_IT615.txt"
    default_average_z_file = "data/average_z_OT800_IT615.txt"

    # Check for help flag
    if '-h' in sys.argv or '--help' in sys.argv:
        print("\nUsage: python compute_pixelmap.py [centroid_file] [geom_file] [average_radius_file] [average_z_file]")
        print("\nOptions:")
        print(f"  centroid_file          Path to the centroid file. Default is {default_centroid_file}")
        print(f"  geom_file              Path to the geometry file. Default is {default_geom_file}")
        print(f"  average_radius_file    Path to the average radius file. Default is {default_average_radius_file}")
        print(f"  average_z_file         Path to the average z file. Default is {default_average_z_file}")
        sys.exit()

    # Determine file paths based on arguments provided
    centroid_file = sys.argv[1] if len(sys.argv) > 1 else default_centroid_file
    geom_file = sys.argv[2] if len(sys.argv) > 2 else default_geom_file
    average_radius_file = sys.argv[3] if len(sys.argv) > 3 else default_average_radius_file
    average_z_file = sys.argv[4] if len(sys.argv) > 4 else default_average_z_file

    det_geom = DetectorGeometry(geom_file, average_radius_file, average_z_file)
    centroid = Centroid(centroid_file)

    # Make output folder if it doesn't exist
    os.makedirs("output/pixelmap/", exist_ok=True)

    printPixelMap(centroid, det_geom)

import os
import sys
import math
import multiprocessing
import pickle

import numpy as np
from tqdm import tqdm

from DetectorGeometry import DetectorGeometry
from Module import Module
from Centroid import Centroid
import LSTMath as sdlmath

ptthresh = 0.8

def get_straight_line_connections_parallel(ref_detid, return_dict={}):
    return_dict[ref_detid] = get_straight_line_connections(ref_detid)
    return True

def get_curved_line_connections_parallel(ref_detid, return_dict={}):
    return_dict[ref_detid] = get_curved_line_connections(ref_detid)
    return True

def get_straight_line_connections(ref_detid):

    # reference module centroid
    centroid = centroidDB.getCentroid(ref_detid)

    # reference module phi position
    refphi = math.atan2(centroid[1], centroid[0])

    # reference module Module instance
    ref_module = Module(ref_detid)

    # reference module layer
    ref_layer = ref_module.layer()

    # reference subdet
    ref_subdet = ref_module.subdet()

    tar_detids_to_be_considered = []
    if ref_subdet == 5:
        tar_detids_to_be_considered += det_geom.getBarrelLayerDetIds(ref_layer + 1)
    else:
        tar_detids_to_be_considered += det_geom.getEndcapLayerDetIds(ref_layer + 1)

    list_of_detids_etaphi_layer_tar = []
    # for tar_detid in tqdm(tar_detids_to_be_considered, desc="looping over target detids"):
    for tar_detid in tar_detids_to_be_considered:
        if sdlmath.module_overlaps_in_eta_phi(
                det_geom.getData()[ref_detid],
                det_geom.getData()[tar_detid],
                refphi,
                0
                ):
            list_of_detids_etaphi_layer_tar.append(tar_detid)
        elif sdlmath.module_overlaps_in_eta_phi(
                det_geom.getData()[ref_detid],
                det_geom.getData()[tar_detid],
                refphi,
                10
                ):
            list_of_detids_etaphi_layer_tar.append(tar_detid)
        elif sdlmath.module_overlaps_in_eta_phi(
                det_geom.getData()[ref_detid],
                det_geom.getData()[tar_detid],
                refphi,
                -10
                ):
            list_of_detids_etaphi_layer_tar.append(tar_detid)

    # Consider barrel to endcap connections if the intersection area is > 0
    if ref_subdet == 5:

        list_of_barrel_endcap_connected_tar_detids = []
        for zshift in [0, 10, -10]:

            ref_polygon = sdlmath.get_etaphi_polygon(det_geom.getData()[ref_detid], refphi, zshift)

            # Check whether there is still significant non-zero area
            for tar_detid in list_of_detids_etaphi_layer_tar:
                tar_polygon = sdlmath.get_etaphi_polygon(det_geom.getData()[tar_detid], refphi, zshift)
                ref_polygon = ref_polygon.difference(tar_polygon)

            # If area is "non-zero" then consider endcap
            tar_detids_to_be_considered = []
            if ref_polygon.area > 0.0001:

                tar_detids_to_be_considered += det_geom.getEndcapLayerDetIds(1)

                # for tar_detid in tqdm(tar_detids_to_be_considered, desc="looping over target detids"):
                for tar_detid in tar_detids_to_be_considered:

                    # If the centroids are far away then don't consider (this was important to exclude incorrect matching when target module is pi away)
                    centroid_target = centroidDB.getCentroid(tar_detid)

                    tarphi = math.atan2(centroid_target[1], centroid_target[0])

                    if abs(sdlmath.Phi_mpi_pi(tarphi - refphi)) > math.pi / 2.:
                        continue

                    tar_polygon = sdlmath.get_etaphi_polygon(det_geom.getData()[tar_detid], refphi, zshift)

                    if ref_polygon.intersects(tar_polygon):
                        list_of_barrel_endcap_connected_tar_detids.append(tar_detid)

        list_of_barrel_endcap_connected_tar_detids = list(set(list_of_barrel_endcap_connected_tar_detids))

        list_of_detids_etaphi_layer_tar += list_of_barrel_endcap_connected_tar_detids

    return list_of_detids_etaphi_layer_tar

def bounds_after_curved(ref_detid, doR=True):
    # Obtaining positive particle helices from centroid and bounds of reference module
    bounds = det_geom.getData()[ref_detid]
    centroid = centroidDB.getCentroid(ref_detid)
    charge = 1
    next_layer_bound_points = []
    theta = math.atan2(math.sqrt(centroid[0]**2 + centroid[1]**2), centroid[2])
    refphi = math.atan2(centroid[1], centroid[0])
    ref_layer = Module(ref_detid).layer()
    ref_subdet = Module(ref_detid).subdet()
    for bound in bounds:
        helix_p10 = sdlmath.construct_helix_from_points(ptthresh, 0, 0,  10, bound[1], bound[2], bound[0], -charge)
        helix_m10 = sdlmath.construct_helix_from_points(ptthresh, 0, 0, -10, bound[1], bound[2], bound[0], -charge)
        helix_p10_pos = sdlmath.construct_helix_from_points(ptthresh, 0, 0,  10, bound[1], bound[2], bound[0], charge)
        helix_m10_pos = sdlmath.construct_helix_from_points(ptthresh, 0, 0, -10, bound[1], bound[2], bound[0], charge)
        # helix_p10 = sdlmath.construct_helix_from_points(1, 0, 0,   0, bound[1], bound[2], bound[0], -charge)
        # helix_m10 = sdlmath.construct_helix_from_points(1, 0, 0,   0, bound[1], bound[2], bound[0], -charge)
        # helix_p10_pos = sdlmath.construct_helix_from_points(1, 0, 0,   0, bound[1], bound[2], bound[0], charge)
        # helix_m10_pos = sdlmath.construct_helix_from_points(1, 0, 0,   0, bound[1], bound[2], bound[0], charge)
        bound_theta = math.atan2(math.sqrt(bound[1]**2 + bound[2]**2), bound[0])
        bound_phi = math.atan2(bound[2], bound[1])

        if ref_subdet == 5:
            tar_layer_radius = det_geom.getBarrelLayerAverageRadius(ref_layer + 1)
            tar_layer_z = det_geom.getEndcapLayerAverageAbsZ(1)

        if ref_subdet == 5:
            if doR:
                tar_layer_radius = det_geom.getBarrelLayerAverageRadius(ref_layer + 1)
                if bound_theta > theta:
                    if sdlmath.Phi_mpi_pi(bound_phi - refphi) > 0:
                        next_point = sdlmath.get_helix_point_from_radius(helix_p10, tar_layer_radius)
                    else:
                        next_point = sdlmath.get_helix_point_from_radius(helix_p10_pos, tar_layer_radius)
                else:
                    if sdlmath.Phi_mpi_pi(bound_phi - refphi) > 0:
                        next_point = sdlmath.get_helix_point_from_radius(helix_m10, tar_layer_radius)
                    else:
                        next_point = sdlmath.get_helix_point_from_radius(helix_m10_pos, tar_layer_radius)
            else:
                tar_layer_z = det_geom.getEndcapLayerAverageAbsZ(1)
                if bound_theta > theta:
                    if sdlmath.Phi_mpi_pi(bound_phi - refphi) > 0:
                        next_point = sdlmath.get_helix_point_from_z(helix_p10, math.copysign(tar_layer_z, helix_p10.lam()))
                    else:
                        next_point = sdlmath.get_helix_point_from_z(helix_p10_pos, math.copysign(tar_layer_z, helix_p10.lam()))
                else:
                    if sdlmath.Phi_mpi_pi(bound_phi - refphi) > 0:
                        next_point = sdlmath.get_helix_point_from_z(helix_m10, math.copysign(tar_layer_z, helix_p10.lam()))
                    else:
                        next_point = sdlmath.get_helix_point_from_z(helix_m10_pos, math.copysign(tar_layer_z, helix_p10.lam()))
        else:
            tar_layer_z = det_geom.getEndcapLayerAverageAbsZ(ref_layer + 1)
            if bound_theta > theta:
                if sdlmath.Phi_mpi_pi(bound_phi - refphi) > 0:
                    next_point = sdlmath.get_helix_point_from_z(helix_p10, math.copysign(tar_layer_z, helix_p10.lam()))
                else:
                    next_point = sdlmath.get_helix_point_from_z(helix_p10_pos, math.copysign(tar_layer_z, helix_p10.lam()))
            else:
                if sdlmath.Phi_mpi_pi(bound_phi - refphi) > 0:
                    next_point = sdlmath.get_helix_point_from_z(helix_m10, math.copysign(tar_layer_z, helix_m10.lam()))
                else:
                    next_point = sdlmath.get_helix_point_from_z(helix_m10_pos, math.copysign(tar_layer_z, helix_m10.lam()))

        next_layer_bound_points.append([next_point[2], next_point[0], next_point[1]])

    return next_layer_bound_points

def get_curved_line_connections(ref_detid):

    # reference module centroid
    centroid = centroidDB.getCentroid(ref_detid)

    # reference module phi position
    refphi = math.atan2(centroid[1], centroid[0])

    # reference module Module instance
    ref_module = Module(ref_detid)

    # reference module layer
    ref_layer = ref_module.layer()

    # reference subdet
    ref_subdet = ref_module.subdet()

    # Target IDs to be considered in the next layer
    tar_detids_to_be_considered = []
    if ref_subdet == 5:
        tar_detids_to_be_considered += det_geom.getBarrelLayerDetIds(ref_layer + 1)
    else:
        tar_detids_to_be_considered += det_geom.getEndcapLayerDetIds(ref_layer + 1)

    next_layer_bound_points = bounds_after_curved(ref_detid)

    list_of_detids_etaphi_layer_tar = []
    # for tar_detid in tqdm(tar_detids_to_be_considered, desc="looping over target detids"):
    for tar_detid in tar_detids_to_be_considered:
        if sdlmath.module_overlaps_in_eta_phi(
                next_layer_bound_points,
                det_geom.getData()[tar_detid],
                refphi,
                0
                ):
            list_of_detids_etaphi_layer_tar.append(tar_detid)
        # elif sdlmath.module_overlaps_in_eta_phi(
        #         next_layer_bound_points,
        #         det_geom.getData()[tar_detid],
        #         refphi,
        #         10
        #         ):
        #     list_of_detids_etaphi_layer_tar.append(tar_detid)
        # elif sdlmath.module_overlaps_in_eta_phi(
        #         next_layer_bound_points,
        #         det_geom.getData()[tar_detid],
        #         refphi,
        #         -10
        #         ):
        #     list_of_detids_etaphi_layer_tar.append(tar_detid)

    # Consider barrel to endcap connections if the intersection area is > 0
    if ref_subdet == 5:

        list_of_barrel_endcap_connected_tar_detids = []
        # for zshift in [0, 10, -10]:
        for zshift in [0]:

            ref_polygon = sdlmath.get_etaphi_polygon(next_layer_bound_points, refphi, zshift)

            # Check whether there is still significant non-zero area
            for tar_detid in list_of_detids_etaphi_layer_tar:
                tar_polygon = sdlmath.get_etaphi_polygon(det_geom.getData()[tar_detid], refphi, zshift)
                ref_polygon = ref_polygon.difference(tar_polygon)

            # If area is "non-zero" then consider endcap
            tar_detids_to_be_considered = []
            if ref_polygon.area > 0.0001:

                tar_detids_to_be_considered += det_geom.getEndcapLayerDetIds(1)

                # for tar_detid in tqdm(tar_detids_to_be_considered, desc="looping over target detids"):
                for tar_detid in tar_detids_to_be_considered:

                    # If the centroids are far away then don't consider (this was important to exclude incorrect matching when target module is pi away)
                    centroid_target = centroidDB.getCentroid(tar_detid)

                    tarphi = math.atan2(centroid_target[1], centroid_target[0])

                    if abs(sdlmath.Phi_mpi_pi(tarphi - refphi)) > math.pi / 2.:
                        continue

                    tar_polygon = sdlmath.get_etaphi_polygon(det_geom.getData()[tar_detid], refphi, zshift)

                    if ref_polygon.intersects(tar_polygon):
                        list_of_barrel_endcap_connected_tar_detids.append(tar_detid)

        list_of_barrel_endcap_connected_tar_detids = list(set(list_of_barrel_endcap_connected_tar_detids))

        list_of_detids_etaphi_layer_tar += list_of_barrel_endcap_connected_tar_detids

    return list_of_detids_etaphi_layer_tar

def write_straight_line_connections(output="output/module_connection_tracing_straight.txt"):
    write_connections(docurved=False,output=output)

def write_curved_line_connections(output="output/module_connection_tracing_curved.txt"):
    write_connections(docurved=True,output=output)

def write_connections(docurved=False, output="output/module_connection_tracing.txt"):
    list_of_detids_etaphi_layer_ref = det_geom.getDetIds(
            lambda x:
            ((Module(x[0]).subdet() == 5 and Module(x[0]).isLower() == 1 and Module(x[0]).layer() != 6) or
            (Module(x[0]).subdet() == 4 and Module(x[0]).isLower() == 1 and Module(x[0]).layer() != 5 and
                not (Module(x[0]).ring() == 15 and Module(x[0]).layer() == 1) and
                not (Module(x[0]).ring() == 15 and Module(x[0]).layer() == 2) and
                not (Module(x[0]).ring() == 12 and Module(x[0]).layer() == 3) and
                not (Module(x[0]).ring() == 12 and Module(x[0]).layer() == 4)
            ))
            # and (Module(x[0]).layer() == 1 and Module(x[0]).rod() == 1 and Module(x[0]).isLower() == 1)
            )
    ref_detid = sorted(list(list_of_detids_etaphi_layer_ref))[0]

    njobs = 32
    manager = multiprocessing.Manager()
    return_dict = manager.dict()
    pool = multiprocessing.Pool(processes=njobs)

    module_map = {}
    for ref_detid in tqdm(list_of_detids_etaphi_layer_ref, desc="looping over ref detids"):
        if docurved:
            job = pool.apply_async(get_curved_line_connections_parallel, args=(ref_detid, return_dict))
            # module_map[ref_detid] = get_curved_line_connections(ref_detid)
        else:
            job = pool.apply_async(get_straight_line_connections_parallel, args=(ref_detid, return_dict))
            # module_map[ref_detid] = get_straight_line_connections(ref_detid)
    pool.close()
    pool.join()

    module_map = dict(return_dict)

    f = open(output, "w")
    print("Writing module connections...")
    for ref_detid in sorted(tqdm(module_map.keys())):
        tar_detids = [str(x) for x in module_map[ref_detid]]
        f.write("{} {} {}\n".format(ref_detid, len(tar_detids), " ".join(tar_detids)))

if __name__ == "__main__":
    # Default file paths
    default_centroid_file = "output/sensor_centroids.txt"
    default_geom_file = "output/sensor_corners.txt"
    default_average_radius_file = "data/average_r_OT800_IT615.txt"
    default_average_z_file = "data/average_z_OT800_IT615.txt"

    # Check for help flag
    if '-h' in sys.argv or '--help' in sys.argv:
        print("\nUsage: python compute_modulemap.py [centroid_file] [geom_file] [average_radius_file] [average_z_file]")
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

    # Setting up detector geometry (centroids and boundaries)
    centroidDB = Centroid(centroid_file)
    dirpath = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    det_geom = DetectorGeometry(geom_file, average_radius_file, average_z_file)
    det_geom.buildByLayer()

    # Make output folder if it doesn't exist
    os.makedirs(os.path.dirname("output/"), exist_ok=True)

    write_curved_line_connections()
    write_straight_line_connections()

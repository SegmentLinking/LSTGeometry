#!/bin/env python

import numpy as np
import os
import math
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
# import mpl_toolkits.mplot3d as a3
import pylab as pl
from DetectorGeometry import DetectorGeometry
from Module import Module
import LSTMath as sdlmath
import LSTDisplay
from Centroid import Centroid
from tqdm import tqdm
import pickle
from matplotlib.collections import LineCollection
import multiprocessing

# ptthresh = 1
ptthresh = 0.8

# Setting up detector geometry (centroids and boundaries)
centroidDB = Centroid("data/centroid.txt")
dirpath = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
det_geom = DetectorGeometry("data/geom.txt", "data/average_radius.txt", "data/average_z.txt")
det_geom.buildByLayer()
sdlDisplay = LSTDisplay.LSTDisplay(det_geom)

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

def write_straight_line_connections(output="data/module_connection_tracing_straight.txt"):
    write_connections(docurved=False,output=output)

def write_curved_line_connections(output="data/module_connection_tracing_curved.txt"):
    write_connections(docurved=True,output=output)

def write_connections(docurved=False, output="data/module_connection_tracing.txt"):

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

def visualize_connection_between_reference_and_target(ref_detid, tar_detids):

    # Retrieve baseline
    # ax_xy = pickle.load(file('/nfs-7/userdata/phchang/detector_layout_matplotlib_pickle/detxy.pickle'))
    # ax_rz = pickle.load(file('/nfs-7/userdata/phchang/detector_layout_matplotlib_pickle/detrz.pickle'))
    ax_xy = pickle.load(file('detxy.pickle'))
    ax_rz = pickle.load(file('detrz.pickle'))

    # Set detector elements
    sdlDisplay.set_detector_xy_collection([ref_detid] + tar_detids)
    sdlDisplay.set_detector_rz_collection([ref_detid] + tar_detids)

    # Draw detector elements
    sdlDisplay.display_detector_xy(ax_xy, color=(0, 0, 0))
    sdlDisplay.display_detector_rz(ax_rz, color=(0, 0, 0))

    # Obtain xy line collections
    ref_centroid = centroidDB.getCentroid(ref_detid)
    tar_centroids = [ centroidDB.getCentroid(x) for x in tar_detids ]

    # List of detids
    detids = []
    detids += tar_detids
    detids.append(ref_detid)

    # Obtain line collection
    segments_rz = []
    segments_xy = []
    centroid_ref = centroidDB.getCentroid(ref_detid)
    ref_x = centroid_ref[0]
    ref_y = centroid_ref[1]
    ref_z = centroid_ref[2]
    ref_rt = math.sqrt(ref_x**2 + ref_y**2)
    for target_detid in tar_detids:
        centroid_target = centroidDB.getCentroid(target_detid)
        target_x = centroid_target[0]
        target_y = centroid_target[1]
        target_z = centroid_target[2]
        target_rt = math.sqrt(target_x**2 + target_y**2)
        segments_rz.append([(ref_z, ref_rt), (target_z, target_rt)])
        segments_xy.append([(ref_x, ref_y), (target_x, target_y)])

    # Create line collection
    lc_rz = LineCollection(segments_rz, colors=(1,0,0), linewidth=0.5, alpha=0.4)
    lc_xy = LineCollection(segments_xy, colors=(1,0,0), linewidth=0.5, alpha=0.4)

    # Draw line collections
    ax_rz.add_collection(lc_rz)
    ax_xy.add_collection(lc_xy)

    # Save the figure
    plt.sca(ax_rz)
    plt.savefig("module_connection_visualization_rz_ref_module_{}.pdf".format(ref_detid))

    plt.sca(ax_xy)
    plt.savefig("module_connection_visualization_xy_ref_module_{}.pdf".format(ref_detid))

    fig, ax = plt.subplots(figsize=(4. * 2,2.*math.pi))

    # reference module centroid
    centroid = np.array(centroidDB.getCentroid(ref_detid))

    # reference module phi position
    refphi = math.atan2(centroid[1], centroid[0])

    ref_polygon = sdlmath.get_etaphi_polygon(
            det_geom.getData()[ref_detid], refphi)

    for tar_detid in tar_detids:
        tar_polygon = sdlmath.get_etaphi_polygon(
                det_geom.getData()[tar_detid],
                refphi
                )

        centroid_target = np.array(centroidDB.getCentroid(tar_detid))

        if np.linalg.norm(centroid_target - centroid) > 100:

            print("printing one connection check case")
            print(refphi)
            print(ref_detid)
            print(det_geom.getData()[ref_detid])
            mod_boundaries = np.array([ sdlmath.get_etaphi([x[1], x[2], x[0]], refphi) for x in det_geom.getData()[ref_detid] ])
            print(mod_boundaries)
            mod_boundaries = np.array([ sdlmath.get_etaphi([x[1], x[2], x[0]], 0) for x in det_geom.getData()[ref_detid] ])
            print(mod_boundaries)
            print(tar_detid)
            print(det_geom.getData()[tar_detid])
            mod_boundaries = np.array([ sdlmath.get_etaphi([x[1], x[2], x[0]], refphi) for x in det_geom.getData()[tar_detid] ])
            print(mod_boundaries)
            mod_boundaries = np.array([ sdlmath.get_etaphi([x[1], x[2], x[0]], 0) for x in det_geom.getData()[tar_detid] ])
            print(mod_boundaries)
            print(ref_polygon.intersects(tar_polygon))
            print(sdlmath.module_overlaps_in_eta_phi(det_geom.getData()[ref_detid], det_geom.getData()[tar_detid], refphi))


    next_layer_bound_points = np.array(bounds_after_curved(ref_detid))
    next_layer_bound_points = np.array([ sdlmath.get_etaphi([x[1], x[2], x[0]], 0) for x in next_layer_bound_points ])
    print(next_layer_bound_points)
    ax.scatter(next_layer_bound_points[0:4,0], next_layer_bound_points[0:4,1])

    sdlDisplay.set_detector_etaphi_collection([ref_detid])
    sdlDisplay.display_detector_etaphi(ax, color=(1,0,0))

    sdlDisplay.set_detector_etaphi_collection(tar_detids)
    sdlDisplay.display_detector_etaphi(ax, color=(0,0,1))
    fig.savefig("module_connection_visualization_etaphi_ref_module_{}.pdf".format(ref_detid))

def visualize_connection(connection_file, ref_detid_to_visualize):

    f = open(connection_file)
    lines = f.readlines()

    sdlDisplay = LSTDisplay.getDefaultLSTDisplay()
    fullLSTDisplay = LSTDisplay.getDefaultLSTDisplay()

    detids = []
    segments_rz = []
    segments_xy = []
    for line in lines:
        ls = line.split()
        ref_detid = int(ls[0])
        nconn = int(ls[1])
        target_detids = []
        if nconn > 0:
            target_detids = [ int(x) for x in ls[2:] ]

        if ref_detid == ref_detid_to_visualize:
            visualize_connection_between_reference_and_target(ref_detid, target_detids)
            break

def visualize_connections(connection_file, ref_detid_to_visualize):

    f = open(connection_file)
    lines = f.readlines()

    sdlDisplay = LSTDisplay.getDefaultLSTDisplay()
    fullLSTDisplay = LSTDisplay.getDefaultLSTDisplay()

    detids = []
    segments_rz = []
    segments_xy = []
    for line in lines:
        ls = line.split()
        ref_detid = int(ls[0])
        nconn = int(ls[1])
        target_detids = []
        if nconn > 0:
            target_detids = [ int(x) for x in ls[2:] ]

        if ref_detid == ref_detid_to_visualize:
            visualize_connection_between_reference_and_target(ref_detid, target_detids)
            break

if __name__ == "__main__":

    write_curved_line_connections()
    write_straight_line_connections()

import sys
import json
import numpy as np
import pandas as pd

def rodrigues_rotation_matrix(axis, theta):
    """
    Calculates the Rodrigues' rotation matrix for rotating a vector around an arbitrary axis.
    
    Parameters:
    axis (numpy.ndarray): The axis around which to rotate (should be a unit vector).
    theta (float): The rotation angle in radians.
    
    Returns:
    numpy.ndarray: The rotation matrix.
    """
    axis = axis / np.linalg.norm(axis)  # Ensure the axis is a unit vector
    K = np.array([
        [0, -axis[2], axis[1]],
        [axis[2], 0, -axis[0]],
        [-axis[1], axis[0], 0]
    ])
    
    R = np.eye(3) + np.sin(theta) * K + (1 - np.cos(theta)) * (K @ K)
    return R

def tangential_rotation_matrix(phi, theta):
    """
    Generates a rotation matrix for rotating around the tangential direction in cylindrical coordinates.
    
    Parameters:
    phi (float): The angular coordinate in cylindrical coordinates (used to define the tangential direction).
    theta (float): The rotation angle around the tangential axis (the "pitching" angle).
    
    Returns:
    numpy.ndarray: The rotation matrix.
    """
    # Tangential unit vector in cylindrical coordinates
    tangential_vector = np.array([-np.sin(phi), np.cos(phi), 0])
    
    # Get the rotation matrix using Rodrigues' rotation formula
    return rodrigues_rotation_matrix(tangential_vector, theta)

def rotation_matrix(tilt_deg, skew_deg, yaw_deg, phi_deg):
    """
    Computes the final rotation matrix based on tilt and phi angles.

    Parameters:
    tilt_deg (float): The tilt angle in degrees.
    skew_deg (float): The skew angle in degrees.
    yaw_deg (float): The yaw angle in degrees.
    phi_deg (float): The phi angle in degrees.
    
    Returns:
    numpy.ndarray: The final rotation matrix.

    Note:
    Only the tilt angles are non-zero for the current geometry. If the other
    angles get used, implement their rotations using the tangential_rotation_matrix
    function above as an example.
    """
    # Check if skew and yaw angles are non-zero
    if skew_deg != 0 or yaw_deg != 0:
        raise ValueError(f"Non-zero skew or yaw angles found. See note in rotation_matrix code."
                         f"Tilt: {tilt_deg}, Skew: {skew_deg}, Yaw: {yaw_deg}, Phi: {phi_deg}")

    # Convert angles from degrees to radians
    tilt_rad = np.radians(tilt_deg)
    phi_rad = np.radians(phi_deg)

    # Rotation around Z-axis that makes the sensor "face towards" the beamline (i.e. towards z-axis)
    # So for example if phi=0 then R is the identity (i.e. already facing), or if phi=90deg
    # then R becomes (x,y,z)->(-y,x,z) so the sensor is rotated 90 degrees to face the beamline
    R_initial = np.array([
        [np.cos(phi_rad), -np.sin(phi_rad), 0],
        [np.sin(phi_rad), np.cos(phi_rad), 0],
        [0, 0, 1]
    ])

    # The tilt angle given in the CSV files is with respect to a module that is facing
    # the beamline, meaning after R_initial is applied. From there we tilt the module according
    # to the rotation below. Note that because this tilt angle is not with respect to the X,Y,Z
    # axes and is instead around an arbitrary axis (defined from the rotation above) we have to apply
    # the Rodrigues' rotation formula
    R_tilt = tangential_rotation_matrix(phi_rad, -tilt_rad)

    R = R_tilt @ R_initial

    return R

def transform_sensor_corners(df_row):
    """
    Calculates the transformed corners of each sensor based on input DataFrame row.

    Parameters:
    df_row (pandas.Series): A row from the DataFrame containing sensor information.
    
    Returns:
    list: A list of transformed corners for the sensor.
    """
    # Extracting the necessary variables from the dataframe row
    det_id = df_row["DetId/U"]
    module_z = df_row[" SensorCenter z(mm)"]
    module_rho = df_row[" SensorCenter rho(mm)"]
    module_phi = df_row[" phi_deg/D"]
    sensor_spacing = df_row[" sensorSpacing_mm/D"]
    sensor_width = df_row[" meanWidth_mm/D"]
    sensor_length = df_row[" length_mm/D"]

    # Convert phi to radians and calculate x, y for the center from rho and phi
    phi_rad = np.radians(module_phi)
    module_x = module_rho * np.cos(phi_rad)
    module_y = module_rho * np.sin(phi_rad)

    # Sensor dimensions used to define corners
    half_width = sensor_width / 2
    half_length = sensor_length / 2
    half_spacing = sensor_spacing / 2

    # Make the module sizes consistent with hit-based method.
    # FIXME: Using the real (smaller) sizes specified by CSV file increases
    # fake rate significantly and lowers efficiency between abs(eta) 1 to 2. 
    width_extension = 50.0 - half_width
    length_extension = (50.0 if half_length > 40 else 25.0) - half_length

    half_width += width_extension
    half_length += length_extension

    # Define the corners of the sensor before rotation
    corners = np.array([
        [-half_spacing, -half_width, -half_length],
        [-half_spacing, -half_width, half_length],
        [-half_spacing, half_width, half_length],
        [-half_spacing, half_width, -half_length],
        [half_spacing, -half_width, -half_length],
        [half_spacing, -half_width, half_length],
        [half_spacing, half_width, half_length],
        [half_spacing, half_width, -half_length]
    ])

    # Apply rotation matrix
    R = rotation_matrix(df_row[" tiltAngle_deg/D"],
                        df_row[" skewAngle_deg/D"],
                        df_row[" yawAngle_deg/D"],
                        df_row[" phi_deg/D"])

    rotated_corners = np.dot(R, corners.T).T
    
    # Apply final coordinate shift
    final_corners = rotated_corners + np.array([module_x, module_y, module_z])

    # Scale the final coordinates to cm from mm
    final_corners_scaled = final_corners / 10

    # Coordinate reorder before saving (x,y,z)->(z,x,y)
    final_corners_adjusted = final_corners_scaled[:, [2, 0, 1]]

    return final_corners_adjusted.tolist()

def assign_corners_to_sensor(module_csv, sensor_csv):
    """
    Assigns each set of four corners to the correct sensor DetID based on the closest centroid.

    Parameters:
    module_csv (pandas.DataFrame): DataFrame containing module information and transformed corners.
    sensor_csv (pandas.DataFrame): DataFrame containing sensor information.
    
    Returns:
    dict: A dictionary with sensor DetID as keys and transformed corners as values.
    """
    transformed_corners_dict = {}

    # Sensor centroid coordinates in Rho, Phi axes
    sensor_rho = sensor_csv[' sensorCenterRho_mm/D']
    sensor_phi = sensor_csv[' phi_deg/D']
    sensor_det = sensor_csv['DetId/i']

    # Sensor centroid coordinates in X,Y,Z axes
    sensor_z = sensor_csv[' sensorCenterZ_mm/D']
    sensor_x = sensor_rho * np.cos(np.radians(sensor_phi))
    sensor_y = sensor_rho * np.sin(np.radians(sensor_phi))

    # Loop through the modules and assign the set of transformed corner coordinates to a sensor
    for index, row in module_csv.iterrows():
        # Module and corresponding sensor DetId's
        module_det_id = row["DetId/U"]
        sensor_det_id_1 = module_det_id + 1
        sensor_det_id_2 = module_det_id + 2

        # Grab transformed corner coordinates and calculate their centroids (mean X,Y,Z position)
        transformed_corners = row['Transformed_Corners']
        centroid_sensor_1 = np.mean(transformed_corners[:4], axis=0)
        centroid_sensor_2 = np.mean(transformed_corners[4:], axis=0)

        # Grab each sensor's true X,Y,Z centroid position from the CSV file
        sensor_index_1 = sensor_csv[sensor_csv['DetId/i'] == sensor_det_id_1].index[0]
        sensor_index_2 = sensor_csv[sensor_csv['DetId/i'] == sensor_det_id_2].index[0]
        sensor_centroid_1 = np.array([sensor_z[sensor_index_1], sensor_x[sensor_index_1], sensor_y[sensor_index_1]])
        sensor_centroid_2 = np.array([sensor_z[sensor_index_2], sensor_x[sensor_index_2], sensor_y[sensor_index_2]])

        # Scale these from mm to cm in order to compare
        sensor_centroid_1 = sensor_centroid_1 / 10
        sensor_centroid_2 = sensor_centroid_2 / 10

        # Calculate distances from true centroids to transformed coordinate centroids and find best match.
        distance_to_sensor_1 = np.linalg.norm(centroid_sensor_1 - sensor_centroid_1)
        distance_to_sensor_2 = np.linalg.norm(centroid_sensor_2 - sensor_centroid_1)
        if distance_to_sensor_1 < distance_to_sensor_2:
            transformed_corners_dict[str(sensor_det_id_1)] = transformed_corners[:4]
            transformed_corners_dict[str(sensor_det_id_2)] = transformed_corners[4:]
        else:
            transformed_corners_dict[str(sensor_det_id_1)] = transformed_corners[4:]
            transformed_corners_dict[str(sensor_det_id_2)] = transformed_corners[:4]
    return transformed_corners_dict

if __name__ == "__main__":
    # Default file paths for module and sensor information
    default_module_info_path = "data/module_info_OT800_IT615.csv"
    default_sensor_info_path = "data/DetId_sensors_list_OT800_IT615.csv"
    default_output_path = "data/sensor_corners.txt"

    # Check for help flag
    if '-h' in sys.argv or '--help' in sys.argv:
        print("\nUsage: python compute_corners_csv.py [module_info_file] [sensor_info_file] [outputfile]")
        print("\nOptions:")
        print(f"  module_info_file  Path to the module information CSV file. Default is {default_module_info_path}")
        print(f"  sensor_info_file  Path to the sensor information CSV file. Default is {default_sensor_info_path}")
        print(f"  outputfile        Path for the output file. Default is {default_output_path}\n")
        sys.exit()

    # Determine input and output file paths based on arguments provided
    module_info_path = sys.argv[1] if len(sys.argv) > 1 else default_module_info_path
    sensor_info_path = sys.argv[2] if len(sys.argv) > 2 else default_sensor_info_path
    output_path = sys.argv[3] if len(sys.argv) > 3 else default_output_path

    # Load the data
    module_csv = pd.read_csv(module_info_path)
    sensor_csv = pd.read_csv(sensor_info_path)

    # Apply the function to each row in the dataframe to get transformed corners
    module_csv['Transformed_Corners'] = module_csv.apply(transform_sensor_corners, axis=1)

    # Get the dictionary with assigned corners
    assigned_corners = assign_corners_to_sensor(module_csv, sensor_csv)

    # Save the dictionary to the specified output file
    with open(output_path, 'w') as f:
        json.dump(assigned_corners, f, indent=4)

    print(f"\nProcessed module info file: {module_info_path}")
    print(f"Processed sensor info file: {sensor_info_path}")
    print(f"Output written to: {output_path}\n")
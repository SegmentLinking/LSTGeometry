## Geometry File Information

The files `DetId_sensors_list.csv` and `module_info.csv` are sourced from Tracker version OT800_IT615:

- URL: [Tracker Version OT800_IT615](https://cms-tklayout.web.cern.ch/cms-tklayout/layouts-work/recent-layouts/OT800_IT615/info.html)

Sources for specific files:
- `allCoordinates.csv` (renamed `module_info.csv`) is available at [OT800_IT615 Layout](https://cms-tklayout.web.cern.ch/cms-tklayout/layouts-work/recent-layouts/OT800_IT615/layout.html)
- `DetId_sensors_list.csv` can be found linked from the homepage of the above URL.
- The `average_r_OT800_IT615.txt` and `average_z_OT800_IT615.txt` files can be taken directly from the table at the top of the [OT800_IT615 Layout](https://cms-tklayout.web.cern.ch/cms-tklayout/layouts-work/recent-layouts/OT800_IT615/layout.html) page. These represent the average r positions of the Barrel layers and the average z positions of the Endcap layers.

## Setting up the relevant environment

    # Download and install Miniconda
    curl -O -L https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b 
    
    # add conda to the end of ~/.bashrc, so relogin after executing this line
    ~/miniconda3/bin/conda init
    
    # (optional) stop conda from activating the base environment on login
    conda config --set auto_activate_base false

    # Add conda-forge as a priority channel for package management
    conda config --add channels conda-forge
    
    # Create a new conda environment with necessary packages
    conda create --name analysisenv uproot pandas matplotlib jupyter graphviz iminuit scipy shapely root
    conda activate analysisenv
    conda install -c plotly plotly=4.14.3
    pip install yahist particle graphviz pydot tqdm

    # Note: After installation, activate your environment with:
    # conda activate analysisenv
    # To deactivate, use:
    # conda deactivate

## Compute Geometry (CSV)

The `compute_geometry.py` file computes both the sensor corner coordinates and centroid coordinates, as well as two orientation files used by the segment linking algorithm, using the `compute_corners.py`, `compute_centroids.py`, and `compute_orientation.py` files respectively. This is the only file that you need to run prior to generating the module maps and pixel maps. Documentation for the individual geometry python files it calls are also given below, but do not need to be run in addition to compute_geometry.

Usage:

    Run: python3 compute_geometry.py for default file paths.

    For custom paths: python3 compute_geometry.py [module_info_file] [sensor_info_file] [outputfile_corners] [outputfile_centroid] [outputfile_tilted_barrel] [outputfile_endcap]
    Default module info file: data/module_info_OT800_IT615.csv
    Default sensor info file: data/DetId_sensors_list_OT800_IT615.csv
    Default output file for corners: output/sensor_corners.txt
    Default output file for centroids: output/sensor_centroids.txt
    Default output file for tilted barrel orientations: output/tilted_barrel_orientation.txt
    Default output file for endcap orientations: output/endcap_orientation.txt

## Computing the Module Maps and Pixel Maps

After running the `compute_geometry.py` file (see above) the following can be run.

Use the scripts:

    python3 compute_pixelmap.py
    python3 compute_modulemap.py

This will place the modulemap output to the output/ directory and the pixelmaps to their own pixelmap directory.

## Convert Outputs to Binary

After creating the module map and pixel map with the above scripts, you can convert the relevant files stored in output/ to a binary format for use in Tracklooper with the following command:

    python3 convert_binary.py

## Compute Centroids (CSV)

The `compute_centroids.py` file computes the centroid coordinates of each sensor using the CSV files in /data

Usage:

    Run: python3 compute_centroids.py for default file paths.

    For custom paths: python3 compute_centroids.py [inputfile] [outputfile]
    Default input: data/DetId_sensors_list_OT800_IT615.csv
    Default output: output/sensor_centroids.txt

Output Format:

    sensor_centroids.txt - [sensor detid], [x coordinate of centroid (cm)], [y coordinate of centorid (cm)], [z coordinate of centroid (cm)], [moduletype (23 (PSP), 24 (PSS), or 25 (TwoS))]

## Compute Corners (CSV)

The `compute_corners.py` file calculates the four corner coordinates of each sensor based on the provided module and sensor CSV files. It uses rotation matrices to account for various rotations of each sensor and outputs the corner coordinates for each sensor.

Usage:

    Run: python3 compute_corners.py for default file paths.

    For custom paths: python3 compute_corners.py [module_info_file] [sensor_info_file] [outputfile]
    Default module info file: data/module_info_OT800_IT615.csv
    Default sensor info file: data/DetId_sensors_list_OT800_IT615.csv
    Default output file: output/sensor_corners.txt

Output Format:

    sensor_corners.txt - "sensor detid": [Z, X, Y coordinates for each of the sensor's four corners (cm)]

## Compute Orientations (CSV)

The `compute_orientation.py` script calculates the orientations (dr/dz and dx/dy slopes) of each relevant sensor based on their corner coordinates. It outputs two files: one for the slopes of tilted barrel sensors and another for the slopes of endcap sensors. Note that only the dxdy slope is given for endcap sensors because dz is always 0 in the current geometry for the endcap sensors. Additionally, for endcap sensors, the centroid phi value also appended to orientation information.

Usage:

    Run: python3 compute_orientation.py for default file paths.

    For custom paths: python3 compute_orientation.py [sensor_corners_file] [centroids_file] [output_tilted_barrel_file] [output_endcap_file]
    Default sensor corners file: output/sensor_corners.txt
    Default centroids file: data/DetId_sensors_list_OT800_IT615.csv
    Default output file for tilted barrel orientations: output/tilted_barrel_orientation.txt
    Default output file for endcap orientations: output/endcap_orientation.txt

Output Format:

    endcap_orientation.txt - [endcap sensor detid] [dx/dy slope of sensor] [centroid phi value of sensor]
    tilted_barrel_orientation.txt - [tilted barrel sensor detid] [dr/dz slope of sensor] [dx/dy slope of sensor]

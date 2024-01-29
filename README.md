## Setting up the relevant environment

    curl -O -L https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b 
    
    # add conda to the end of ~/.bashrc, so relogin after executing this line
    ~/miniconda3/bin/conda init
    
    # stop conda from activating the base environment on login
    conda config --set auto_activate_base false
    conda config --add channels conda-forge
    
    conda create --name analysisenv uproot pandas matplotlib jupyter graphviz iminuit scipy shapely root
    conda activate analysisenv # To be repeated when re-doing this study
    conda install -c plotly plotly=4.14.3
    pip install yahist particle graphviz pydot tqdm
    
    # Optional steps
    git clone git@github.com:SegmentLinking/LSTStudy.git
    cd LSTStudy
    
    jupyter notebook --no-browser
    
    # remember the port number and do something like the following on your local computer
    ssh -N -f -L localhost:$port:localhost:$port uaf-$uaf.t2.ucsd.edu

## Geometry File Information

The files `DetId_sensors_list.csv` and `module_info.csv` are sourced from Tracker version OT800_IT615:

- URL: [Tracker Version OT800_IT615](https://cms-tklayout.web.cern.ch/cms-tklayout/layouts-work/recent-layouts/OT800_IT615/info.html)

Sources for specific files:
- `allCoordinates.csv` (renamed `module_info.csv`) is available at [OT800_IT615 Layout](https://cms-tklayout.web.cern.ch/cms-tklayout/layouts-work/recent-layouts/OT800_IT615/layout.html)
- `DetId_sensors_list.csv` can be found linked from the homepage of the above URL.
- The `average_r_OT800_IT615.txt` and `average_z_OT800_IT615.txt` files can be taken directly from the table at the top of the [OT800_IT615 Layout](https://cms-tklayout.web.cern.ch/cms-tklayout/layouts-work/recent-layouts/OT800_IT615/layout.html) page. These represent the average r positions of the Barrel layers and the average z positions of the Endcap layers.

## Compute Geometry (CSV)

The `compute_geometry.py` file computes both the sensor corner coordinates and centroid coordinates using the `compute_corners.py` and `compute_centroids.py` files respectively.

Usage:

    Run `python3 compute_geometry.py` for default file paths.

    For custom paths: `python3 compute_geometry.py [module_info_file] [sensor_info_file] [outputfile_corners] [outputfile_centroid]`
    Replace `[module_info_file]`, `[sensor_info_file]`, `[outputfile_corners]`, and `[outputfile_centroid]` with your specific file paths.
    Default module info file: `data/module_info_OT800_IT615.csv`
    Default sensor info file: `data/DetId_sensors_list_OT800_IT615.csv`
    Default output file for corners: `output/sensor_corners.txt`
    Default output file for centroids: `output/sensor_centroids.txt`

## Compute Centroids (CSV)

The `compute_centroids.py` file computes the centroid coordinates of each sensor using the CSV files in /data

Usage:

    Run `python3 compute_centroids.py` for default file paths.

    For custom paths: `python3 compute_centroids.py [inputfile] [outputfile]`
    Replace `[inputfile]` and `[outputfile]` with your specific file paths.
    Default input: `data/DetId_sensors_list_OT800_IT615.csv`
    Default output: `output/sensor_centroids.txt`

## Compute Corners (CSV)

The `compute_corners.py` file calculates the four corner coordinates of each sensor based on the provided module and sensor CSV files. It uses rotation matrices to account for various rotations of each sensor and outputs the corner coordinates for each sensor.

Usage:

    Run `python3 compute_corners.py` for default file paths.

    For custom paths: `python3 compute_corners.py [module_info_file] [sensor_info_file] [outputfile]`
    Replace `[module_info_file]`, `[sensor_info_file]`, and `[outputfile]` with your specific file paths.
    Default module info file: `data/module_info_OT800_IT615.csv`
    Default sensor info file: `data/DetId_sensors_list_OT800_IT615.csv`
    Default output file: `output/sensor_corners.txt`

## Computing the module maps and pixel maps

After running the `compute_geometry.py` file (see above) the following can be run.

Use the script:

    compute_pixelmap.py
    compute_modulemap.py

This will place the modulemap output to the output/ directory and the pixelmaps to their own pixelmap directory.
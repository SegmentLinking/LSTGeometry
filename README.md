## Creating geometry related files

Take a PU200 ttbar sample so that we have enough hits to perform some simple fits to gather data on the geometry.


## Conda for jupyter notebook

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

## Printing Hits

First step of generating the module maps and pixel maps.
Saves the hits in the outer tracker container to a file.
One can use the tracking ntuple from the CMSSW output.

    print_hits.py

    Usage:

    Run the script with the following command:

    python print_hits.py [inputfile] [outputfile] [nevents]

    inputfile: Optional. The path to the ROOT file. If not specified, a default path is used.
    outputfile: Optional. The path to the output text file. If not specified, defaults to hits.txt.
    nevents: Optional. The maximum number of events to process. If not specified, defaults to 100.

    The script writes the hit data in the following format:

    x: [x-coordinate] y: [y-coordinate] z: [z-coordinate] detId: [detector ID] moduleType: [module type]

    For example:

    x: 109.938 y: -7.6565 z: 101.17 detId: 443363422 moduleType: 25
    x: 110.055 y: -6.24833 z: 101.17 detId: 443363422 moduleType: 25
    x: 109.921 y: -7.86728 z: 106.194 detId: 443363422 moduleType: 25
    x: 109.656 y: -5.61791 z: 115.163 detId: 443363425 moduleType: 25
    x: 109.517 y: -7.29516 z: 115.163 detId: 443363425 moduleType: 25
    x: 109.481 y: -5.54475 z: 115.163 detId: 443363426 moduleType: 25
    x: 109.342 y: -7.22648 z: 115.163 detId: 443363426 moduleType: 25
    ...

## Compute Centroids (CSV)

Alternative way of computing the sensor centroids using the CSV files in /data/ directly 
instead of using the hit information computed in the previous step.

Usage:

    Run `python compute_centroid_csv.py` for default file paths.

    For custom paths: `python compute_centroid_csv.py [inputfile] [outputfile]`
    Replace `[inputfile]` and `[outputfile]` with your specific file paths.
    Default input: `data/DetId_sensors_list_OT800_IT615.csv`
    Default output: `data/centroid.txt`

## Computing the Centroids of the modules

Using the printed hits in the above format in a txt file, point the file path in the `compute_centroid.py` script and run it.
It will then write an output to `data/centroid.txt`.

## Computing the endcap/tilted module orientation

Use the scripts. Change the input path to the correct txt file.

    fit_endcap_module_orientations.py
    fit_tilted_module_orientations.py

They will place output to:

    data/endcap_orientation.txt
    data/tilted_orientation.txt

## Computing the module 4 corner bounds

Use the script:

    compute_geometry.py

This will place output to:

    data/geom.txt

## Computing the module maps between the module to another module in a different layer and with the pixel seeds

Use the script:

    compute_pixelmap.py
    compute_connection.py

If one looks at the ```if __name__ == "__main__":``` area of the script, there are two main functions to be called:

This will place output to:

    data/module_connection_tracing.txt

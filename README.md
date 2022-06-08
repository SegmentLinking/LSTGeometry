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

## Printing Hits

The script contains some example of how one wants to print out the hits in the outer tracker container.
One can use the tracking ntuple from the CMSSW output.

    print_hit.C

They will be printed out in the following format. Pipe them to a txt file.

    x: 28.0882 y: 1.38098 z: -129.049 detId: 411309061
    x: 27.0796 y: 5.87351 z: -129.049 detId: 411309061
    x: 27.5126 y: 1.26449 z: -129.049 detId: 411309061
    x: 27.3692 y: 1.23166 z: -129.049 detId: 411309061
    x: 26.6629 y: 5.69108 z: -129.049 detId: 411309061
    x: 27.1555 y: 1.6433 z: -129.049 detId: 411309061
    x: 26.7019 y: 4.50759 z: -129.049 detId: 411309061
    x: 27.2155 y: 0.326476 z: -129.049 detId: 411309061
    x: 26.6088 y: 1.34409 z: -129.049 detId: 411309061
    x: 26.1426 y: 4.2874 z: -129.049 detId: 411309061
    ...

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

## Computing the module maps between the module to another module in a different layer

Use the script:

    compute_connection.py

If one looks at the ```if __name__ == "__main__":``` area of the script, there are two main functions to be called:

This will place output to:

    data/module_connection_tracing.txt

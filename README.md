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
    
    conda create --name analysisenv uproot pandas matplotlib jupyter graphviz iminuit scipy shapely
    conda activate analysisenv
    conda install -c plotly plotly=4.14.3
    pip install yahist particle graphviz pydot tqdm
    
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
    x: 26.4285 y: 0.606807 z: -129.049 detId: 411309061
    x: 26.1897 y: 1.17647 z: -129.049 detId: 411309061
    x: 25.9786 y: 2.50985 z: -129.049 detId: 411309061
    x: 25.2019 y: 7.41372 z: -129.049 detId: 411309061
    x: 25.4848 y: 4.68945 z: -129.049 detId: 411309061
    x: 25.0562 y: 7.39571 z: -129.049 detId: 411309061
    x: 24.7423 y: 8.43947 z: -129.049 detId: 411309061
    x: 25.8026 y: -0.130182 z: -129.049 detId: 411309061
    x: 25.7744 y: 0.0476022 z: -129.049 detId: 411309061
    x: 25.2652 y: 3.26253 z: -129.049 detId: 411309061
    x: 25.2457 y: 3.38599 z: -129.049 detId: 411309061
    x: 25.0063 y: 4.89715 z: -129.049 detId: 411309061
    x: 25.082 y: 3.48156 z: -129.049 detId: 411309061
    x: 24.8974 y: 4.64703 z: -129.049 detId: 411309061
    x: 24.5079 y: 7.10638 z: -129.049 detId: 411309061
    x: 25.0075 y: 3.01415 z: -129.049 detId: 411309061
    x: 24.373 y: 6.08267 z: -129.049 detId: 411309061
    x: 24.3174 y: 6.4333 z: -129.049 detId: 411309061
    x: 24.031 y: 7.30421 z: -129.049 detId: 411309061
    x: 24.6573 y: 2.41196 z: -129.049 detId: 411309061
    x: 24.6151 y: 2.67863 z: -129.049 detId: 411309061
    x: 23.6772 y: 8.59982 z: -129.049 detId: 411309061
    x: 24.7705 y: 0.75932 z: -129.049 detId: 411309061
    x: 23.6473 y: 7.85092 z: -129.049 detId: 411309061
    x: 24.4369 y: 0.989975 z: -129.049 detId: 411309061
    x: 24.2469 y: 2.19002 z: -129.049 detId: 411309061
    x: 23.9645 y: 3.97279 z: -129.049 detId: 411309061
    x: 23.9614 y: 3.99255 z: -129.049 detId: 411309061
    x: 24.3124 y: 0.838627 z: -129.049 detId: 411309061
    x: 23.33 y: 7.04131 z: -129.049 detId: 411309061

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

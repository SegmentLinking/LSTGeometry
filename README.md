# Creating geometry related files

Take a PU200 ttbar sample so that we have enough hits to perform some simple fits to gather data on the geometry.


# Conda for jupyter notebook

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
    pip install yahist particle graphviz pydot
    
    git clone git@github.com:SegmentLinking/LSTStudy.git
    cd LSTStudy
    
    jupyter notebook --no-browser
    
    # remember the port number and do something like the following on your local computer
    ssh -N -f -L localhost:$port:localhost:$port uaf-$uaf.t2.ucsd.edu
# LSTGeometry

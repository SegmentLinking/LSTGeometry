## (Deprecated!) Creating geometry related files from hits

Take a PU200 ttbar sample so that we have enough hits to perform some simple fits to gather data on the geometry.

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

## Computing the Centroids of the modules

Using the printed hits in the above format in a txt file, point the file path in the `compute_centroids.py` script and run it.
It will then write an output to `output/sensor_centroids.txt`.

## Computing the endcap/tilted module orientation

Use the scripts. Change the input path to the correct txt file.

    fit_endcap_module_orientations.py
    fit_tilted_module_orientations.py

They will place output to:

    output/endcap_orientation.txt
    output/tilted_orientation.txt

## Computing the sensor 4 corner bounds

Use the script:

    compute_geometry.py

This will place output to:

    output/sensor_corners.txt

## Computing the average r (Barrel) and z (Endcap) values of each layer

Use the script:

    compute_average_layer_geom.py

The above script stores the average radius for each layer in the barrel and the average z value for each layer in the endcap

This will place the output to:

    output/average_radius.txt
    output/average_z.txt

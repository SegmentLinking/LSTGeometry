import struct
import os
import glob

def ensure_dir(directory):
    """Ensure that the output directory exists."""
    if not os.path.exists(directory):
        os.makedirs(directory)

def convert_variable_columns_to_binary(input_file, output_file):
    """Convert files with variable number of columns to binary format."""
    with open(input_file, 'r') as f, open(output_file, 'wb') as out:
        for line in f:
            parts = line.split()
            data = [int(parts[0]), int(parts[1])] + [int(x) for x in parts[2:]]
            fmt_str = 'II' + 'I' * data[1]
            out.write(struct.pack(fmt_str, *data))

def convert_fixed_columns_to_binary(input_file, output_file, fmt_str):
    """Convert files with fixed number of columns to binary format."""
    with open(input_file, 'r') as f, open(output_file, 'wb') as out:
        for line in f:
            parts = line.replace(',', ' ').split()
            if fmt_str == 'IfffI':  # Specific case for sensor_centroids.txt
                data = [int(parts[0]), float(parts[1]), float(parts[2]), float(parts[3]), int(parts[4])]
            else:  # General case for other files
                data = [int(parts[0])] + [float(x) for x in parts[1:]]
            out.write(struct.pack(fmt_str, *data))

def convert_pixelmap_files(input_dir, output_dir):
    """Convert all files in the pixelmap directory."""
    pixelmap_input_dir = os.path.join(input_dir, 'pixelmap')
    pixelmap_output_dir = os.path.join(output_dir, 'pixelmap')
    ensure_dir(pixelmap_output_dir)  # Create the pixelmap output directory
    files = glob.glob(os.path.join(pixelmap_input_dir, '*.txt'))
    for file in files:
        output_file = os.path.join(pixelmap_output_dir, os.path.basename(file).replace('.txt', '.bin'))
        convert_variable_columns_to_binary(file, output_file)

# Input + Output directory
input_dir = 'output/'
output_dir = 'bin_output/'

# Create the output directory
ensure_dir(output_dir)

# Convert module map to binary
convert_variable_columns_to_binary(input_dir+'module_connection_tracing_merged.txt', output_dir+'module_connection_tracing_merged.bin')

# Convert all files in the pixelmap directory
convert_pixelmap_files(input_dir, output_dir)

# Convert fixed column files to binary
convert_fixed_columns_to_binary(input_dir+'endcap_orientation.txt', output_dir+'endcap_orientation.bin', 'Iff')
convert_fixed_columns_to_binary(input_dir+'sensor_centroids.txt', output_dir+'sensor_centroids.bin', 'IfffI')
convert_fixed_columns_to_binary(input_dir+'tilted_barrel_orientation.txt', output_dir+'tilted_barrel_orientation.bin', 'Iff')
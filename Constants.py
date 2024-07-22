# Constants for Barrel and Endcap subdet type
BARREL = 5
ENDCAP = 4

# Constants for module types
PSP = 23
PSS = 24
TwoS = 25

# Default value for slopes if undefined (taken from Tracklooper Repo)
DEFAULT_SLOPE = float("inf")

# Minimum pT threshold for the module map and pixel map
PTTHRESH = 0.8

# Magnetic field constant (Tesla)
B = 3.8112

# Constant used to get curvature of track (GeV/c/(T*cm))
k2Rinv1GeVf = 0.00299792458

# Default file names
DEFAULT_MODULE_INFO_FILE = "data/module_info_OT800_IT615.csv"
DEFAULT_SENSOR_INFO_FILE = "data/DetId_sensors_list_OT800_IT615.csv"
DEFAULT_AVERAGE_R_FILE = "data/average_r_OT800_IT615.txt"
DEFAULT_AVERAGE_Z_FILE = "data/average_z_OT800_IT615.txt"

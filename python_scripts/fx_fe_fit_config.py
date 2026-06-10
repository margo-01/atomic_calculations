# fit_config.py

DATA_DIR = '/home/margo/core_valence_fe_scripts/python_scripts/data'  # Directory where the .out files are located

BASE_DIR = '/home/margo/core_valence_fe_scripts/python_scripts/data/output_data/all_e'  # Base directory where the .xy files are located
OUT_CSV = "GL_test.csv" # Name of the output csv file where the fitted parameters and metrics will be saved

# Cromer–Mann parameterisation  model is preferred for fx, no constrain on parameters because valence fx has negative regions, set NG to 5 and NC to 1
# 
# Gaussian/Lorentzian model is preferred for the parametrisation of fe, the a,b,c,d, parameters are constrained to be positive, set NL and NG to 3, and NC to 0


# Number of terms
NL = 3   # Lorentzian-like terms, only used for GL
NG = 3   # Gaussian terms
NC = 0 # Number of constant terms in the model
Nparams = 2 * NL + 2 * NG + NC


# Optimisation settings
N_STARTS = 2
MAX_NFEV = 10000

# Bounds for all parameters
LOWER_BOUND = 0.0001
UPPER_BOUND = 10.0

# Elements to fit
ELEMENTS = ["C", "N", "O", "Si", "P"]






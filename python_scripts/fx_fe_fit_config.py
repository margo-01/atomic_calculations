# fit_config.py

DATA_DIR = '/home/margo/core_valence_fe_scripts/python_scripts/data'  # Directory where the .out files are located

#BASE_DIR = '/home/margo/core_valence_fe_scripts/python_scripts/data/output_data/fx_all_e'  # Base directory where the .xy files are located
#OUT_CSV = "M_test.csv" # Name of the output csv file where the fitted parameters and metrics will be saved

# Cromer–Mann parameterisation  model is preferred for fx, no constrain on parameters because valence fx has negative regions, set NG to 5 and NC to 1
# 
# Gaussian/Lorentzian model is preferred for the parametrisation of fe, the a,b,c,d, parameters are constrained to be positive, set NL and NG to 3, and NC to 0
#model_type = "M"  # "G" for Gaussian model, "GL" for Gaussian/Lorentzian model, "M" for a mixed unconstrained model 
#data_type = "all_e"  # "core", "val", or "all_e" for all-electron data
# Number of terms
#NL = 3   # Lorentzian-like terms, only used for GL
#NG = 3   # Gaussian terms
#NC = 0 # Number of constant terms in the model
#Nparams = 2 * NL + 2 * NG + NC


# Optimisation settings
N_STARTS = 1000
MAX_NFEV = 10000

# Bounds for all parameters
#LOWER_BOUND = 0.0001
#UPPER_BOUND = 10.0

#FX = True  # Set to True if fitting fx data, False for fe data
# Elements to fit
ELEMENTS = ["H","He","Li","Be","B","C","N","O","F","Ne", "Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca", "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe"]


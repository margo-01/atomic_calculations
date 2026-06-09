# fit_config.py

BASE_DIR = "/data/margo/atomic_calculations/fe_elements/RSCAN"
OUT_CSV = "GL_params_all_elements_Be_core.csv"

# Cromer–Mann parameterisation  model is preferred for fx, a>0 , b>0, c>0, set NG to 5 and NC to 1
# 
# Gaussian/Lorentzian model is preferred for the parametrisation of fe, the a,b,c,d, parameters are constrained to be positive, set NL and NG to 3, and NC to 0


# Number of terms
NL = 3   # Lorentzian-like terms, only used for GL
NG = 3   # Gaussian terms
NC = 0 # Number of constant terms in the model
Nparams = 2 * NL + 2 * NG + NC


# Optimisation settings
N_STARTS = 500
MAX_NFEV = 10000

# Bounds for all parameters
LOWER_BOUND = 0.0001
UPPER_BOUND = 10.0

# Elements to fit
ELEMENTS = ["Be"]






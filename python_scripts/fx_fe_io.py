import numpy as np
from pathlib import Path

#the loading functions rely on the same file naming and structure as in the github repo

def load_element_data(data_dir, element):
    # this function is used to import all-e core and valence data for a selected element
    #either fx or electron density 
    data_dir = Path(data_dir)
    element = element
    recip_grid = np.linspace(0, 52.917720860738790, 2001) * 1.889 / (2 * np.pi)
    all_e_data = np.loadtxt(data_dir/'all_e'/f"{element}_*.out")
    core_data = np.loadtxt(data_dir/'core'/f"{element}_*.out")
    val_data = np.loadtxt(data_dir/'val'/f"{element}_*.out")
    
    all_e = np.array([recip_grid, all_e_data])
    core = np.array([recip_grid, core_data])
    val = np.array([recip_grid, val_data])
   
    return all_e, core, val


def load_fx_norm():
    pass

def load_fe_data():
    pass

def load_fe_norm():
    pass

def save_scattering_data():
    pass

def save_plot(fig, filename, dpi=300):
    fig.savefig(filename, dpi=dpi, bbox_inches="tight")
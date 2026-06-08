import numpy as np
from pathlib import Path

#the loading functions rely on the same file naming and structure as in the github repo

def load_fx_data(data_dir, element):
    # this function is used to import all-e core and valence data for a selected element
    data_dir = Path(data_dir)
    element = element

    recip_grid = np.linspace(0, 52.917720860738790, 2001) * 1.889 / (2 * np.pi)


    all_e_data = np.loadtxt(data_dir/'all_e'/f"{element}_recip_iam_radial_charge.out")
    core_data = np.loadtxt(data_dir/'core'/f"{element}_recip_core_radial_charge.out")
    val_data = np.loadtxt(data_dir/'val'/f"{element}_recip_val_radial_charge.out")

    all_e = np.array([recip_grid, all_e_data])
    core = np.array([recip_grid, core_data])
    val = np.array([recip_grid, val_data])
    
    return all_e, core, val

def load_rho_data(data_dir, element):
    # this function is used to import all-e core and valence data for a selected element
    data_dir = Path(data_dir)
    element = element
    all_e_data = np.loadtxt(data_dir/'all_e'/f"{element}_rho_ae.out")
    core_data = np.loadtxt(data_dir/'core'/f"{element}_rho_core.out")
    val_data = np.loadtxt(data_dir/'val'/f"{element}_rho_val.out")


   
    return all_e_data, core_data, val_data


def load_fx_norm():
    pass

def load_fe_data():
    pass

def load_fe_norm():
    pass

def save_fe_data(out_dir, element, f_e, kind = None):

    f_e_out = np.column_stack((f_e))
    if kind is not None:
        filename = f"{element}_fe_{kind}.xy"
    else:
        filename = f"{element}_fe.xy"
    np.savetxt(out_dir/filename, f_e_out)


def save_plot(fig, filename, element=  None,dpi=300):
    filename = Path(filename)
    if element is not None:
        filename = filename.with_name(f"{filename.stem}_{element}{filename.suffix or '.png'}")
    fig.savefig(filename, dpi=dpi, bbox_inches="tight")

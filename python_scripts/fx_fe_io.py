import numpy as np
from pathlib import Path
import pandas as pd
import csv
import os


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

def load_f_e_x_data(filename):
    data = np.loadtxt(filename)

    x = data[:, 0]
    y = data[:, 1]

    return x, y


def load_fe_norm():
    pass

def save_fe_data(out_dir, element, recip_grid, f_e, kind=None):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    f_e_out = np.column_stack((recip_grid, f_e))
    if kind is not None:
        filename = f"{element}_fe_{kind}.xy"
    else:
        filename = f"{element}_fe.xy"
    np.savetxt(out_dir/filename, f_e_out)

def save_fx_data(out_dir, element, recip_grid, f_x, kind=None):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    f_x_out = np.column_stack((recip_grid, f_x))
    if kind is not None:
        filename = f"{element}_fx_{kind}.xy"
    else:
        filename = f"{element}_fx.xy"
    np.savetxt(out_dir/filename, f_x_out)


def save_plot(fig, filename, element=  None,dpi=300):
    filename = Path(filename)
    filename.parent.mkdir(parents=True, exist_ok=True)
    if element is not None:
        filename = filename.with_name(f"{filename.stem}_{element}{filename.suffix or '.png'}")
    fig.savefig(filename, dpi=dpi, bbox_inches="tight")



def GL_write_param_csv(BASE_DIR, OUT_CSV, row, NL, NG):
    out_path = os.path.join(BASE_DIR, OUT_CSV)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    header = ["type", "symbol", "Z", "RMSE", "chi2_red"] + \
            [f"a{i+1}" for i in range(NL)] + [f"b{i+1}" for i in range(NL)] + \
            [f"c{i+1}" for i in range(NG)] + [f"d{i+1}" for i in range(NG)]

    if not os.path.exists(out_path):
        with open(out_path, "w", newline="") as f:
            csv.writer(f).writerow(header)   
    with open(out_path, "a", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(row)


def G_write_param_csv(BASE_DIR, OUT_CSV, row, NG, NC):
    out_path = os.path.join(BASE_DIR, OUT_CSV)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    header = ["type", "symbol", "Z", "RMSE", "chi2_red"] + \
            [f"a{i+1}" for i in range(NG)] + [f"b{i+1}" for i in range(NG)] + \
            [f"c{i+1}" for i in range(NC)] 

    if not os.path.exists(out_path):
        with open(out_path, "w", newline="") as f:
            csv.writer(f).writerow(header)   
    with open(out_path, "a", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(row)

def M_write_param_csv(BASE_DIR, OUT_CSV, row, NL, NG, NC):
    out_path = os.path.join(BASE_DIR, OUT_CSV)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    header = ["type", "symbol", "Z", "RMSE", "chi2_red"] + \
            [f"a{i+1}" for i in range(NL)] + [f"b{i+1}" for i in range(NL)] + \
            [f"c{i+1}" for i in range(NG)] + [f"d{i+1}" for i in range(NG)] + [f"cc{i+1}" for i in range(NC)] 

    if not os.path.exists(out_path):
        with open(out_path, "w", newline="") as f:
            csv.writer(f).writerow(header)   
    with open(out_path, "a", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(row)

def load_GL_param(data_dir, element):
    df = pd.read_csv(data_dir)#
    coeffs_by_Z = {}
    
    for element, group in df.groupby(element):
        n_el = group['n_el'].to_numpy(dtype=float)
        a_vals = group[['a1', 'a2', 'a3', 'a4']].to_numpy(dtype=float)
        b_vals = group[['b1', 'b2', 'b3', 'b4']].to_numpy(dtype=float)
        c_vals = group[['c1', 'c2', 'c3', 'c4']].to_numpy(dtype=float)
        d_vals = group[['d1', 'd2', 'd3', 'd4']].to_numpy(dtype=float)
        cc_vals = group[['cc1']].to_numpy(dtype=float)

        coeffs_by_Z[element] = {'n_el': n_el, 'a': a_vals, 'b': b_vals, 'c': c_vals, 'd': d_vals, 'cc1': cc_vals}

    a_vals = np.asarray(coeffs_by_Z['a'], dtype=float)
    b_vals = np.asarray(coeffs_by_Z['b'], dtype=float)
    c_vals = np.asarray(coeffs_by_Z['c'], dtype=float)
    d_vals = np.asarray(coeffs_by_Z['d'], dtype=float)
    cc_vals = np.asarray(coeffs_by_Z['cc1'], dtype=float)

    return a_vals, b_vals, c_vals, d_vals, cc_vals


def load_G_param(data_dir, Z):
    pass













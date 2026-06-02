import os
import glob
import numpy as np
from scipy.optimize import least_squares
import csv
import matplotlib.pyplot as plt

BASE_DIR = "/data/margo/atomic_calculations/fe_elements/RSCAN_norm_val"   
OUT_CSV = "GL_params_all_elements_fx_val_norm.csv"

 
NL = 3
NG = 3
Nparams = 2 * NL + 2 * NG
n_starts = 500
max_nfev = 10000

lb = [0.0001] * Nparams
ub = [10.0] * Nparams

header = ["type", "symbol", "Z", "RMSE", "chi2_red"] + \
         [f"a{i+1}" for i in range(NL)] + [f"b{i+1}" for i in range(NL)] + \
         [f"c{i+1}" for i in range(NG)] + [f"d{i+1}" for i in range(NG)]

if not os.path.exists(os.path.join(BASE_DIR, OUT_CSV)):
    with open(os.path.join(BASE_DIR, OUT_CSV), "w", newline="") as f:
        csv.writer(f).writerow(header)   

elements = [ "H","He","Li","Be","B","C","N","O","F","Ne", "Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca", "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe"]#,"Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"]



def model(p, q, NL=NL, NG=NG):
    a = p[0:NL]
    b = p[NL:2*NL]
    c = p[2*NL:2*NL+NG]
    d = p[2*NL+NG:2*NL+2*NG]
    # shapes: a[:,None] / (q**2 + b[:,None]) -> (NL, len(q)), sum over axis=0 -> (len(q),)
    y_param = (a[:, None] / (q**2 + b[:, None])).sum(axis=0) + (c[:, None] * np.exp(-d[:, None] * q**2)).sum(axis=0)
    return y_param

def residuals(p, x, y):
    return model(p, x) - y

def random_start(Nparams, lb, ub):
    lb = np.array(lb); ub = np.array(ub)
    r = np.random.rand(Nparams)

    return lb + r * (ub - lb)

def fit_one(x, y, lb, ub, n_starts=200, max_nfev=2000):
    best = None
    best_cost = np.inf
    # make arrays for bounds
    lb_arr = np.array(lb)
    ub_arr = np.array(ub)
    for _ in range(n_starts):
        p0 = random_start(Nparams, lb_arr, ub_arr)
        try:
            res = least_squares(residuals, p0, args=(x, y), method='lm',
                                xtol=1e-9, ftol=1e-9, gtol=1e-9, max_nfev=max_nfev)
        except Exception as e:
            # sometimes LS may fail; skip that start
            # print("least_squares exception:", e)
            continue
        if res.cost < best_cost:
            best_cost = res.cost
            best = res
    return best

for i in elements:
    folder_name = os.path.join(BASE_DIR, i)
    os.chdir(folder_name)
    data_val = np.loadtxt('fx_val_norm.xy')
    
    x = data_val[:, 0]
    y_val = data_val[:,1]

    
    best_val = fit_one(x, y_val, lb, ub, n_starts=n_starts, max_nfev=max_nfev)


    params_val = best_val.x
    fitted_val = model(params_val, x)
    resid_val = y_val - fitted_val

    rmse_val = np.sqrt(np.mean(resid_val**2))


    chi2_red_val = np.sum(resid_val**2) / max(1, len(y_val) - Nparams )

    

    row_val = ["valence", i , elements.index(i) + 1, rmse_val, chi2_red_val] + [float(v) for v in params_val]

    with open(os.path.join(BASE_DIR, OUT_CSV), "a", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(row_val )


    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 6), sharex=True,
                               gridspec_kw={'height_ratios': [3, 1]})

    ax1.plot(x, y_val, 'o', ms=4, label='DFT data')
    ax1.plot(x, fitted_val, '-', label='Best fit')
    ax1.set_ylabel("Scattering factor")
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, 4)

    ax2.plot(x, resid_val, 'o', ms=3)
    ax2.axhline(0, color='k', lw=1)
    ax2.set_xlabel("x")
    ax2.set_ylabel("Residual")
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0,4)

    plt.tight_layout()
    plt.savefig("fx_fit_res_val_norm.png")
    plt.close(fig)
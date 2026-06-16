import numpy as np
from scipy.optimize import least_squares
from pathlib import Path
import warnings

import fx_fe_fit_config as config
import fx_fe_io as io
import fx_fe_plots as plots

def model_mixed(p, q, NL=0, NG=0, NC=0):
    idx = 0
    y_param = np.zeros_like(q, dtype=float)

    if NL > 0:
        a = p[idx:idx + NL]
        idx += NL

        b = p[idx:idx + NL]
        idx += NL

        y_param += (
            a[:, None]/ (q**2 + b[:, None])
        ).sum(axis=0)

    # Gaussian terms
    if NG > 0:
        c = p[idx:idx + NG]
        idx += NG

        d = p[idx:idx + NG]
        idx += NG

        y_param += (
            c[:, None] * np.exp(-d[:, None] * q**2)
        ).sum(axis=0)

    # Constant term
    if NC > 0:
        e = p[idx:idx + NC]
        idx += NC

        y_param += e.sum()

    return y_param

def modelGL(p, q, NL=0, NG=0):
    a = p[0:NL]
    b = p[NL:2*NL]
    c = p[2*NL:2*NL+NG]
    d = p[2*NL+NG:2*NL+2*NG]
    
    y_param = (a[:, None]**2 / (q**2 + b[:, None]**2)).sum(axis=0) + (c[:, None]**2 * np.exp(-d[:, None]**2 * q**2)).sum(axis=0)
    return y_param

def modelG(p, q, NG=0, NC=0):
    a = p[0:NG]
    b = p[NG:2*NG]
    const = p[2*NG:2*NG+NC]
    y_param = (a[:, None]**2 * np.exp(-b[:, None]**2 * q**2)).sum(axis=0) + const.sum()
    return y_param

def residuals(model, p, x, y):
    return model(p, x) - y

def random_start(Nparams, lb, ub):
    lb = np.array(lb); ub = np.array(ub)
    r = np.random.rand(Nparams)

    return lb + r * (ub - lb)

def fit_data(model, x, y, Nparams, lb, ub, n_starts=200, max_nfev=2000):
    best = None
    best_cost = np.inf
    # make arrays for bounds
    lb_arr = np.array(lb)
    ub_arr = np.array(ub)
    for _ in range(n_starts):
        p0 = random_start(Nparams, lb_arr, ub_arr)
        try:
            res = least_squares(lambda p: residuals(model, p, x, y), p0,
                                method='lm', xtol=1e-9, ftol=1e-9, gtol=1e-9,
                                max_nfev=max_nfev)
        except Exception as e:
            # sometimes LS may fail; skip that start
            # print("least_squares exception:", e)
            continue
        if res.cost < best_cost:
            best_cost = res.cost
            best = res
    return best


def gof_parameters(residuals, y, Nparams):
    chi2_red = np.sum(residuals**2) / max(1, len(y) - Nparams )
    rmse = np.sqrt(np.mean(residuals**2))

    return rmse, chi2_red


def parametrisation_cfg(): # data type can be all_e, core, or val 
    data_type = config.data_type  # Change this to "core" or "val" as needed
    BASE_DIR, OUT_CSV = config.BASE_DIR, config.OUT_CSV
    BASE_DIR = Path(BASE_DIR)
    NL, NG, NC = config.NL, config.NG, config.NC

    Nparams = config.Nparams

    LOWER_BOUND = config.LOWER_BOUND
    UPPER_BOUND = config.UPPER_BOUND

    ELEMENTS = config.ELEMENTS
    MODEL_TYPE = config.model_type
    FX = config.FX

    Z_elements = {
        "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8,
        "F": 9, "Ne": 10, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15,
        "S": 16, "Cl": 17, "Ar": 18, "K": 19, "Ca": 20, "Sc": 21, "Ti": 22,
        "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29,
        "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36,
        "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43,
        "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
        "Sb": 51, "Te": 52, "I": 53, "Xe": 54,
    }
    if MODEL_TYPE == "G":
        model = lambda p, q: modelG(p, q, NG=NG, NC=NC)
        if NL != 0:
            warnings.warn("NL is set to a non-zero value, but the Gaussian model does not use Lorentzian terms. NL will be ignored.")
    elif MODEL_TYPE == "GL":
        model = lambda p, q: modelGL(p, q, NL=NL, NG=NG)
        if NC != 0:
            warnings.warn("NC is set to a non-zero value, but the Gaussian/Lorentzian model does not use constant terms. NC will be ignored.")
    elif MODEL_TYPE == "M":
        model = lambda p, q: model_mixed(p, q, NL=NL, NG=NG, NC=NC)
        warnings.warn("Using mixed model with NL, NG, and NC terms. Coefficients are not constrained to be positive. Ensure that the parameters are set correctly for your data type.")


    for element in ELEMENTS:
        Z = Z_elements.get(element)

        if FX:
            x , y = io.load_f_e_x_data(filename=BASE_DIR /  f"{element}_fx_{data_type}.xy") 
        else:
            x , y = io.load_f_e_x_data(filename=BASE_DIR /  f"{element}_fe_{data_type}.xy")

        best_fit = fit_data(model, x, y, Nparams, lb=[LOWER_BOUND]*Nparams, ub=[UPPER_BOUND]*Nparams,
                            n_starts=config.N_STARTS, max_nfev=config.MAX_NFEV)
        params = best_fit.x
        fitted_final = model(params, x)
        residuals = fitted_final - y 
        rmse, chi2_red = gof_parameters(residuals, y, Nparams=Nparams)


# checking the model type to determine which csv writing function to use, and how to format the data row accordingly, then saving the data
        if MODEL_TYPE == "G":
            data_row = [data_type, element, Z, rmse, chi2_red] + [float(v) for v in params**2]
            io.G_write_param_csv(BASE_DIR/ "Fit_Tables", OUT_CSV, data_row, NG, NC)
        elif MODEL_TYPE == "GL":
            data_row = [data_type, element, Z, rmse, chi2_red] + [float(v) for v in params**2] # squaring the parameters to ensure positivity, as the model depends on the squares of the parameters
            io.GL_write_param_csv(BASE_DIR/ "Fit_Tables", OUT_CSV, data_row, NL, NG)
        elif MODEL_TYPE == "M":
            data_row = [data_type, element, Z, rmse, chi2_red] + [float(v) for v in params] # squaring the parameters to ensure positivity, as the model depends on the squares of the parameters
            io.M_write_param_csv(BASE_DIR/ "Fit_Tables", OUT_CSV, data_row, NL, NG, NC)
        else:
            warnings.warn("Invalid model type. Please choose 'G' for Gaussian model or 'GL' for Gaussian/Lorentzian model.")

        
        if FX:
            io.save_plot(plots.fe_x_param_plot(x, y, fitted_final, residuals), BASE_DIR / "Fit_plots" / "fx_fit.png", element, dpi=300)
        else:
            io.save_plot(plots.fe_x_param_plot(x, y, fitted_final, residuals), BASE_DIR / "Fit_plots" / "fe_fit.png", element, dpi=300)


        print(  f"Finished fitting {data_type} data for element {element}. RMSE: {rmse:.4f}, chi2_red: {chi2_red:.4f}")



def parametrisation(data_type, BASE_DIR, OUT_CSV, NL, NG, NC, LOWER_BOUND, UPPER_BOUND, MODEL_TYPE, FX): # data type can be all_e, core, or val 

    BASE_DIR = Path(BASE_DIR)

    Nparams = 2 * NL + 2 * NG + NC


    ELEMENTS = config.ELEMENTS


    Z_elements = {
        "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8,
        "F": 9, "Ne": 10, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15,
        "S": 16, "Cl": 17, "Ar": 18, "K": 19, "Ca": 20, "Sc": 21, "Ti": 22,
        "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29,
        "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36,
        "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43,
        "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
        "Sb": 51, "Te": 52, "I": 53, "Xe": 54,
    }
    if MODEL_TYPE == "G":
        model = lambda p, q: modelG(p, q, NG=NG, NC=NC)
        if NL != 0:
            warnings.warn("NL is set to a non-zero value, but the Gaussian model does not use Lorentzian terms. NL will be ignored.")
    elif MODEL_TYPE == "GL":
        model = lambda p, q: modelGL(p, q, NL=NL, NG=NG)
        if NC != 0:
            warnings.warn("NC is set to a non-zero value, but the Gaussian/Lorentzian model does not use constant terms. NC will be ignored.")
    elif MODEL_TYPE == "M":
        model = lambda p, q: model_mixed(p, q, NL=NL, NG=NG, NC=NC)
        warnings.warn("Using mixed model with NL, NG, and NC terms. Coefficients are not constrained to be positive. Ensure that the parameters are set correctly for your data type.")


    for element in ELEMENTS:
        Z = Z_elements.get(element)

        if FX:
            x , y = io.load_f_e_x_data(filename=BASE_DIR /  f"{element}_fx_{data_type}.xy") 
        else:
            x , y = io.load_f_e_x_data(filename=BASE_DIR /  f"{element}_fe_{data_type}.xy")

        best_fit = fit_data(model, x, y, Nparams, lb=[LOWER_BOUND]*Nparams, ub=[UPPER_BOUND]*Nparams,
                            n_starts=config.N_STARTS, max_nfev=config.MAX_NFEV)
        if best_fit is None:
            print(f"Skipping {element}: no successful fit")
            continue

        params = best_fit.x
        fitted_final = model(params, x)
        residuals = fitted_final - y 
        rmse, chi2_red = gof_parameters(residuals, y, Nparams=Nparams)


# checking the model type to determine which csv writing function to use, and how to format the data row accordingly, then saving the data
        if MODEL_TYPE == "G":
            data_row = [data_type, element, Z, rmse, chi2_red] + [float(v**2) for v in params[:NG]] + [float(v**2) for v in params[NG:2*NG]] \
    + [float(v) for v in params[2*NG:2*NG+NC]] # avoid squaring the constant, it is not constrained to be positive
            io.G_write_param_csv(BASE_DIR/ "Fit_Tables", OUT_CSV, data_row, NG, NC)
        elif MODEL_TYPE == "GL":
            data_row = [data_type, element, Z, rmse, chi2_red] + [float(v) for v in params**2] # squaring the parameters to ensure positivity, as the model depends on the squares of the parameters
            io.GL_write_param_csv(BASE_DIR/ "Fit_Tables", OUT_CSV, data_row, NL, NG)
        elif MODEL_TYPE == "M":
            data_row = [data_type, element, Z, rmse, chi2_red] + [float(v) for v in params] # 
            io.M_write_param_csv(BASE_DIR/ "Fit_Tables", OUT_CSV, data_row, NL, NG, NC)
        else:
            warnings.warn("Invalid model type. Please choose 'G' for Gaussian model or 'GL' for Gaussian/Lorentzian model.")

        io.save_plot(plots.fe_x_param_plot(x, y, fitted_final, residuals), BASE_DIR / "Fit_plots" / "fe_fit.png", element, dpi=300)

        print(  f"Finished fitting {data_type} data for element {element}. RMSE: {rmse:.4f}, chi2_red: {chi2_red:.4f}")







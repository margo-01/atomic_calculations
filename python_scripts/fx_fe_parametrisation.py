import numpy as np
from scipy.optimize import least_squares

import fx_fe_fit_config as config

def model(p, q, NL=0, NG=0, NC=0):
    idx = 0
    y_param = np.zeros_like(q, dtype=float)

    if NL > 0:
        a = p[idx:idx + NL]
        idx += NL

        b = p[idx:idx + NL]
        idx += NL

        y_param += (
            a[:, None]**2 / (q**2 + b[:, None]**2)
        ).sum(axis=0)

    # Gaussian terms
    if NG > 0:
        c = p[idx:idx + NG]
        idx += NG

        d = p[idx:idx + NG]
        idx += NG

        y_param += (
            c[:, None]**2 * np.exp(-d[:, None]**2 * q**2)
        ).sum(axis=0)

    # Constant term
    if NC > 0:
        e = p[idx:idx + NC]
        idx += NC

        y_param += e.sum()

    return y_param


def residuals(p, x, y):
    return model(p, x) - y

def random_start(Nparams, lb, ub):
    lb = np.array(lb); ub = np.array(ub)
    r = np.random.rand(Nparams)

    return lb + r * (ub - lb)

def fit_data(x, y, lb, ub, n_starts=200, max_nfev=2000):
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


def gof_parameters(residuals, y):
    chi2_red = np.sum(residuals**2) / max(1, len(y) - Nparams )
    rmse = np.sqrt(np.mean(residuals**2))

    return rmse, chi2_red

    pass

def parametrisation():
    pass








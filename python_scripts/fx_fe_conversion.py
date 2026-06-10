# Conversion of x-ray form factors to electron scattering factors
# Constants and conventions adopted from Lobato et al. https://doi.org/10.1107/S205327331401643X
import numpy as np
from pathlib import Path
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import fx_fe_io as io 
import fx_fe_plots as plots
import fx_fe_fit_config as config

def mott_bethe(x, f_x, Z):

    a_0 = 0.529177210544
    x = np.asarray(x, dtype=float)
    f_x = np.asarray(f_x, dtype=float)
    f_e = (1.0 / (2.0 * np.pi**2 * a_0)) * ((Z - f_x) / (x**2))
    
    return f_e

def convert_fx_to_fe(data_dir, element):
    main_dir = Path(data_dir)
    output_dir = main_dir / "output_data"
    output_dir.mkdir(parents=True, exist_ok=True)
    for subdir in ("all_e", "core", "val"):
        (output_dir / subdir).mkdir(parents=True, exist_ok=True)

    output_dir_plots = main_dir / "plots"
    output_dir_plots.mkdir(parents=True, exist_ok=True)
    a_0 = 0.529177210544

    rho_ae, rho_core, rho_val = io.load_rho_data(main_dir/'rho_data/', element)
    fx_all_e, fx_core, fx_val = io.load_fx_data(main_dir/'scattering_data', element)

    recip_grid = fx_all_e[0]

    Z = round(np.trapz(rho_ae[:, 1], rho_ae[:, 0]))
    Z_core = round(np.trapz(rho_core[:, 1], rho_core[:, 0]))
    Z_val = round(np.trapz(rho_val[:, 1], rho_val[:, 0]))

    print('Element: {}, Z: {}, Z_core: {}, Z_val: {}'.format(element, Z, Z_core, Z_val))

    fig = plots.fx_ae_core_val_plot(recip_grid, fx_core[1], fx_val[1], fx_all_e[1], Z)
    io.save_plot(fig, output_dir_plots / "fx_ae_core_val_plot.png", element=element)
    plt.close(fig)



    r2_mean_ae = np.trapz(rho_ae[:, 0]**2 * rho_ae[:, 1], rho_ae[:, 0]) / np.trapz(rho_ae[:, 1], rho_ae[:, 0])
    r2_mean_core = np.trapz(rho_core[:, 0]**2 * rho_core[:, 1], rho_core[:, 0]) / np.trapz(rho_core[:, 1], rho_core[:, 0])
    r2_mean_val = np.trapz(rho_val[:, 0]**2 * rho_val[:, 1], rho_val[:, 0]) / np.trapz(rho_val[:, 1], rho_val[:, 0])

# the following lines are needed to ensure the correct behavior of fe at g=0, where the Mott-Bethe formula diverges. 
# The value of fe at g=0 is determined by the mean square radius of the electron density distribution, as shown in Lobato et al. https://doi.org/10.1107/S205327331401643X

    fe_ae = mott_bethe(recip_grid, fx_all_e[1], Z)
    fe_core = mott_bethe(recip_grid, fx_core[1], Z_core)
    fe_val = mott_bethe(recip_grid, fx_val[1], Z_val)

    fe_ae[0] = (1 / 3) * a_0 * Z * r2_mean_ae
    fe_core[0] = (1 / 3) * a_0 * Z_core * r2_mean_core
    fe_val[0] = (1 / 3) * a_0 * Z_val * r2_mean_val

# however in practice the analytically calculated fe values at g=0 don't match too well so interpolation is used 

    f1 = interp1d(recip_grid, fe_ae, kind="cubic", fill_value="extrapolate")
    fe_ae[0] = f1(recip_grid[0])

    f2 = interp1d(recip_grid, fe_core, kind="cubic", fill_value="extrapolate")
    fe_core[0] = f2(recip_grid[0])

    f3 = interp1d(recip_grid, fe_val, kind="cubic", fill_value="extrapolate")
    fe_val[0] = f3(recip_grid[0])

    fig = plots.fe_ae_core_val_plot(recip_grid, fe_core, fe_val, fe_ae)
    io.save_plot(fig, output_dir_plots / "fe_ae_core_val_plot.png", element=element)
    plt.close(fig)
    
    #saving fe data for later use
    io.save_fe_data(output_dir/'all_e', element, recip_grid, fe_ae, kind="all_e")
    io.save_fe_data(output_dir/'core', element, recip_grid, fe_core, kind="core")
    io.save_fe_data(output_dir/'val', element, recip_grid, fe_val, kind="val")

def mb_conversion():
    ELEMENTS = config.ELEMENTS
    DATA_DIR = config.DATA_DIR
    for element in ELEMENTS:
        convert_fx_to_fe(DATA_DIR, element)    





























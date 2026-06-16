import numpy as np

import fx_fe_conversion as conversion
import fx_fe_parametrisation as param



conversion.mb_conversion(normalised=True)




data_type = ["all_e", "core", "val", "val_norm"]
data_type_fx = ["all_e", "val_norm"]


base_dir_fe = ["/home/margo/core_valence_fe_scripts/python_scripts/data/output_data/fe_all_e", "/home/margo/core_valence_fe_scripts/python_scripts/data/output_data/fe_core", "/home/margo/core_valence_fe_scripts/python_scripts/data/output_data/fe_val", "/home/margo/core_valence_fe_scripts/python_scripts/data/output_data/fe_val_norm"] #"/home/margo/core_valence_fe_scripts/python_scripts/data/output_data/fe_all_e"

base_dir_fx = ["/home/margo/core_valence_fe_scripts/python_scripts/data/output_data/fx_all_e", "/home/margo/core_valence_fe_scripts/python_scripts/data/output_data/fx_val_norm"]

out_csv_fe = ["GL_fe_all_e.csv", "GL_fe_core.csv", "GL_fe_val.csv", "GL_fe_val_norm.csv"]

out_csv_fx = ["GL_fx_all_e.csv", "GL_fx_val_norm.csv" ]


for i in range(len(data_type)):
    param.parametrisation(data_type[i], base_dir_fe[i], out_csv_fe[i], 3, 3, 0, 0.0001, 10.00, "GL", FX=False)



for i in range(len(data_type_fx)):
    param.parametrisation(data_type_fx[i], base_dir_fx[i], out_csv_fx[i], 4, 4, 1, 0.0001, 10.00, "M", FX=True)









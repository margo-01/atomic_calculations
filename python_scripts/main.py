import numpy as np

import fx_fe_conversion as conversion

data_dir = '/home/margo/core_valence_fe_scripts/python_scripts/data'

elements = ["C", "N", "O", "Si", "P"]


for element in elements:
    conversion.convert_fx_to_fe(data_dir, element)





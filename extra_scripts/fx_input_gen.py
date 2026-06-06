import os 

elements = ['K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe']#['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']

def write_input(element):
    filenames = [f"{element}.cell", f"{element}.param" ]
    lattice = [(5.0, 0.0, 0.0),(0.0, 5.0, 0.0),(0.0, 0.0, 5.0),]
    
    species = [f"""{name}"""]
    positions_frac = [(0.5, 0.5, 0.5)]
    kgrid = (1, 1, 1)
    
    with open(filenames[0], "w") as f:
        f.write("%BLOCK lattice_cart\n")
        f.write("   ANG\n")
        for v in lattice:
            f.write(f"       {v[0]:.14f}        {v[1]:.14f}        {v[2]:.14f}\n")
        f.write("%ENDBLOCK lattice_cart\n\n")

        f.write("%BLOCK positions_frac\n")
        for s, p in zip(species, positions_frac):
            f.write(f"   {s}                {p[0]:.15f}       {p[1]:.15f}       {p[2]:.15f}\n")
        f.write("%ENDBLOCK positions_frac\n\n")

        f.write(f"kpoint_mp_grid :   {kgrid[0]}  {kgrid[1]}   {kgrid[2]}\n\n")
        f.write("SYMMETRY_GENERATE\n\n")
        f.write("FIX_ALL_CELL : TRUE\n")
    
    with open(filenames[1], "w") as f1:
        f1.write("cut_off_energy = 100 eV \n")
        f1.write("xc_functional : RSCAN \n")
        f1.write("write_checkpoint : none \n")
        f1.write("%block devel_code\n")
        f1.write("CALCULATE_XRD_SF:\n")
        f1.write(":ENDCALCULATE_XRD_SF\n")
        f1.write("%endblock devel_code\n")

        
base_dir = '/data/margo/atomic_calculations/fe_elements/RSCAN'

for name in elements:
    folderpath = os.path.join(base_dir, name)
    os.makedirs(folderpath)
    #now writing castep input files into the relevant directory
    os.chdir(folderpath)
    write_input(name)
    os.chdir(base_dir)
    
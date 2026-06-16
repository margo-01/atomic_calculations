# atomic_calculations



supported parametrisation models:

$f_{e,\mathrm{param}}(q) = \sum^{N_G} a_i\cdot e^{-\left(\frac{x-b_i}{c_i}\right)^2} + cc$ (a,b > 0)

$f_{e,\mathrm{param}}(q) = \sum^{N_L}\frac{a_i}{q^2+b_i} + \sum^{N_G}c_i e^{-d_iq^2}$ (a,b,c,d > 0)

$f_{e,\mathrm{param}}(q) = \sum^{N_L}\frac{a_i}{q^2+b_i} + \sum^{N_G}c_i e^{-d_iq^2} + cc$ (a,b,c,d, cc not constrained)

in the /data directory:
- atomic electron density on radial logarithmic grid
- reciprocal atomic electron density on radial grid (atomic x-ray form factor)
- .cvs table(s) with parametrisation coefficients


scattering vector defined as: $g = 0.5 \cdot s = \frac{\sin(\theta)}{2\cdot \lambda}$(see Lobato et al. 10.1107/S205327331401643X )

THIS DIFFERS FROM THE CONVENTIONAL INTERNATIONSL TABLES DEFINITION! (to compare the Gaussian model with tabulated values from the international tables - divide the b parameter in the exponent by 4)






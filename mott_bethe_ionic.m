function [f_e] = mott_bethe(x, f_x, Z, dZ)
a_0 = 0.529177210544;
f_e = 0.02393366*(((Z-f_x)./(x.^2)) + dZ./(x.^2));

end
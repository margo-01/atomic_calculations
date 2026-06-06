function [f_e] = mott_bethe(x, f_x)
a_0 = 0.529177210544;
f_e = (1/(2*pi^2*a_0))*((1-f_x)./(x.^2));
end
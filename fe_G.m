function f = fe_G(x,a,b)

x = x(:)';     % 1 × N
a = a(:);      % M × 1
b = b(:);      % M × 1

f = sum(a .* exp(-b .* (x.^2)/4), 1); %WATCH OUT HERE: the parameters in many references are defined for s, not g as in Lobato (s = 2g), hence division by 4

end
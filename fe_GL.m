function f = fe_GL(x, a, b, c, d)
    a = a(:);
    b = b(:);
    c = c(:);
    d = d(:);
    f = sum(c .* exp(-d .* x.^2) + a ./ (x.^2 + b) );
end
function f = fe_lobato(x, a, b)
    a = a(:);
    b = b(:);
    f = a.' *((2+b.*x.^2)./((1+b.*x.^2)).^2)
end
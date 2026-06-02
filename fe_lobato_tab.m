function f = fe_lobato_tab(x, Z)
    coeffs = readtable("/data/margo/atomic_calculations/Lobato_coeffs_dataframe.csv");
    a = table2array(coeffs(Z, 4:8))';
    b = table2array(coeffs(Z, 9:13));
    a = a(:);
    b = b(:);
    f = a.' *((2+b.*x.^2)./((1+b.*x.^2)).^2)
    
end
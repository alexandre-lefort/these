function c = symtbx_get_sympoly_coeffs(poly)

%% Give coefficients of a symbolic polynomial

    syms s;
    [coefs,orders] = coeffs(poly,s);

    syms aux;
    orders2 = subs(poly2sym(orders,aux),{aux},{1});
    orders3 = sym2poly(orders2);
    d = length(orders3) - 1;
    orders4 =sym(orders3); 
    idx = 1;
    for ii=1:d+1
        if (orders3(ii) == 1)
            orders4(ii) = coefs(idx);
            idx = idx +1;
        end
    end
	
    c = orders4;

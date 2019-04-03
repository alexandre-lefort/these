

%% Compute Bm coefficients for Krasovskii H2 Criteria


function coef = compute_bm_eq(poly_a, poly_b, p, k)

	a        = symtbx_poly_coeffs(poly_a,p);
	da       = symtbx_poly_degree(poly_a,p);

	b        = symtbx_poly_coeffs(poly_b,p);
	db       = symtbx_poly_degree(poly_b,p);
	
	coef     = simplify((b(k)/b(1) - a(k)/a(1))^2);
	
	idx      = 1;
	
    while(((k + idx) <= db+1) && ((k - idx) >= 1))
        sgn   = 2*mod(idx,2)-1;
		
		cplus = simplify((b(k + idx)/b(1) - a(k + idx)/a(1)));
		cmins = simplify((b(k - idx)/b(1) - a(k - idx)/a(1)));
		
        coef  = simplify(coef - sgn * 2 *cplus * cmins);
        idx   = idx + 1;
	end
	
	coef = simplify(b(1)^2*coef);
	
end
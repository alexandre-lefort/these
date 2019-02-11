

%% Compute Bm coefficients for Krasovskii H2 Criteria


function coef = compute_bm(poly, p, m)

	c        = symtbx_poly_coeffs(poly,p);
	d        = symtbx_poly_degree(poly,p);
    coef     = c(m)^2;
	idx      = 1;
	
    while(((m + idx) <= d+1) && ((m - idx) >= 1))
        sgn   = 2*mod(idx,2)-1;
        coef  = coef - sgn * 2 * c(m + idx) * c(m - idx);
        idx   = idx + 1;
	end
	
end
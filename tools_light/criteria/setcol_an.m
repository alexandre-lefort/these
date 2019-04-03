

%% Build Hurwitz matrix (auxiliar function) : second method


function M = setcol_an(M, poly, p, col)

    poly_deg   = symtbx_poly_degree(poly,p);
	poly_coeff = symtbx_poly_coeffs(poly,p);
	
    M(1,col) = poly_coeff(poly_deg);
    M(2,col) = poly_coeff(poly_deg + 1);
	
    if (poly_deg > 2)
        for ii = 3:size(M,1)
            M(ii,col) = 0;
		end
	end
	
end
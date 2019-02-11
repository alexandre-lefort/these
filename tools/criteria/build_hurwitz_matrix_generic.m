

%% Build Hurwitz matrix : first method (Popov : 41.14)


function M = build_hurwitz_matrix_old(poly, p)

    d = symtbx_poly_degree(poly,p);
	c = symtbx_poly_coeffs(poly,p);
	
    M = sym(zeros(d));
	
    for jj=1:d
	
	    even_idx = 2*(jj);
		
        for ii=1:d
			
			idx = (even_idx - ii);
			
            if (idx > d)
                M(ii,jj) = 0;
			elseif (idx >= 0)
				 M(ii,jj) = c(idx+1);
			end

		end
	end
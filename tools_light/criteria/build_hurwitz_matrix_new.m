

%% Build Hurwitz matrix : second method


function M = build_hurwitz_matrix_new(poly, p)

    d = symtbx_poly_degree(poly,p);
	c = symtbx_poly_coeffs(poly,p);
	
    M = sym(zeros(d));
	
    for jj=1:d
	
	    even_idx = 2*(jj-1)+1;
		
        for ii=1:d
		
			sgn = -(2*mod(ii+jj,2)-1);
			
			idx = (d - even_idx + ii);
			
            if (idx > d)
                M(ii,jj) = 0;
			elseif (idx >= 0)
				 M(ii,jj) = sgn*simplify(c(idx+1));
			end

		end
	end

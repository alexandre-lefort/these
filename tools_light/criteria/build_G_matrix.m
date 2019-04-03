

%% Build G matrix : auxiliar for impulsion criteria


function G = build_G_matrix(M, poly, p)

    G  = M;
	
    d  = symtbx_poly_degree(poly,p);
	c  = symtbx_poly_coeffs(poly,p);

    d
        
    n = size(M,2); 


    for ii=1:n
        if (2*ii <= d+2)
		    G(1,ii) = simplify(c(2*(ii)-1));
        else
            G(1,ii) = 0;
        end
    end

    if (d==2)
        G(1,3) = G(1,2);
        G(1,2) = G(1,1);
        G(1,1) = 0;
        G(1,:) = -G(1,:);
    end

end
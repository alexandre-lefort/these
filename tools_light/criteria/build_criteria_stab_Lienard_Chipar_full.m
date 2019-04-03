%% Give stability criteria with Lienard Chipar expression :
%% https://en.wikipedia.org/wiki/Li%C3%A9nard%E2%80%93Chipart_criterion


function criteria_coefs = build_criteria_stab_Lienard_Chipar_full(den, p)
	
    den_deg   = symtbx_poly_degree(den,p);
	den_coeff = symtbx_poly_coeffs(den,p);

    M_hurwitz = build_hurwitz_matrix_old(den,p);
	criteria_coefs = {};

    for ii = 1:den_deg+1
        criteria_coefs{end+1} = den_coeff(ii);
    end
	
    for ii = 1:den_deg
        criteria_coefs{end+1} = det(symtbx_minor_matrix(M_hurwitz,ii));
    end

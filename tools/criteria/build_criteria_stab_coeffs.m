function criteria_coefs = build_criteria_stab_coeffs(den)

%% Give stability criteria with Lienard Chipar expression :
%% https://en.wikipedia.org/wiki/Li%C3%A9nard%E2%80%93Chipart_criterion
	
    den_deg   = symtbx_get_sympoly_degree(den);
    den_coeff = symtbx_get_sympoly_coeffs(den);

    criteria_coefs = {};

    for ii = 1:length(den_coeff)
        criteria_coefs{ii} = den_coeff(ii);
    end

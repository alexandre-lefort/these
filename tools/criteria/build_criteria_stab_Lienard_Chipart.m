%% Give stability criteria with Lienard Chipar expression :
%% https://en.wikipedia.org/wiki/Li%C3%A9nard%E2%80%93Chipart_criterion


function criteria_coefs = build_criteria_stab_Lienard_Chipart(den)

    den_deg   = symtbx_get_sympoly_degree(den);
    den_coeff = simplify(symtbx_get_sympoly_coeffs(den));

    M_hurwitz = build_hurwitz_matrix_old(den);
    criteria_coefs = {};

    mode = mod(den_deg,2);

    if (mode == 0)
        ii_c = (den_deg+1):-2:1;
        ii_m = (den_deg)  :-2:1;
    else
        ii_c = (den_deg+1):-2:1;
        ii_m = (den_deg-1):-2:1;
    end

    s_c = length(ii_c);
    s_m = length(ii_m);

    for ii = 1:s_c
        idx = ii_c(ii);
        %% criteria_coefs{end+1} = symtbx_horner(den_coeff(idx),'nodegree');
        criteria_coefs{end+1} = den_coeff(idx);
    end
	
    for ii = 1:s_m
        idx = ii_m(ii);
        minor = symtbx_minor_matrix(M_hurwitz,idx);
        %% criteria_coefs{end+1} = symtbx_horner(det(minor),'nodegree');
        criteria_coefs{end+1} = det(minor);
    end

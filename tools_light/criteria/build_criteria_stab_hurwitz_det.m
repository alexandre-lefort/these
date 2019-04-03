

%% Give stability criteria with Routh Hurwitz criteria :
%% https://en.wikipedia.org/wiki/Routh%E2%80%93Hurwitz_stability_criterion


function [criteria_coefs] = build_criteria_stab_hurwitz_det(det_stab, p)
	
    den_deg   = symtbx_poly_degree(det_stab,p);

    disp(strcat('degree = ', num2str(den_deg) ));

    %% Normalize

    det_stab = simplify(det_stab);
    
    M_hurwitz = build_hurwitz_matrix_old(det_stab,p);
    M_hurwitz
    
	criteria_coefs = cell(1,den_deg);

    for ii = 1:den_deg
        m = symtbx_minor_matrix(M_hurwitz,ii);
        criteria_coefs{ii} = simplify(det(m));
    end


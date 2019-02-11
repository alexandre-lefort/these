

%% Compute Analytic H2 impulsion norm for the fpoly system


function [h2_imp] = build_criteria_h2_imp(fpoly, p)

	[num_poly ,den_poly]  = numden(fpoly);
	num_deg   = symtbx_poly_degree(num_poly,p);
	den_deg   = symtbx_poly_degree(den_poly,p);

	if (num_deg >= den_deg)
		error('degree too big');
	end

	num_poly_min = subs(num_poly,{p},{-p});

	poly_prod = expand(num_poly*num_poly_min);

    M_hurwitz = simplify(build_hurwitz_matrix_old(den_poly,p));
    G_matrix  = simplify(build_G_matrix(M_hurwitz,poly_prod,p));

	num_deg   = symtbx_poly_degree(num_poly,p);
	den_coeff = symtbx_poly_coeffs(den_poly,p);

    sgn    = (-1)^num_deg;
	det_G  = simplify(det(G_matrix));

	det_H  = simplify(det(M_hurwitz));
	den_C  = simplify(den_coeff(1));

    h2_imp = simplify(sgn*det_G/(det_H*2*den_C));

	%% digitsOld = digits(5);
	%% h2_imp = vpa(h2_imp,4);
	%% h2_imp = simplify(h2_imp);

	[num, den] = numden(h2_imp);

	h2_imp = recurse_horner_new(num)/recurse_horner_new(den);

end
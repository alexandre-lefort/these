

%% Compute Analytic H2 step norm for the fpoly system


function [h2_step] = build_criteria_h2_step(fpoly,p)

	[poly_b,poly_a] = numden(fpoly);
	
    M_hurwitz = build_hurwitz_matrix_new(poly_a,p);
		
    db = symtbx_poly_degree(poly_b,p);
    da = symtbx_poly_degree(poly_a,p);

	b = simplify(symtbx_poly_coeffs(poly_b,p));	
	a = simplify(symtbx_poly_coeffs(poly_a,p));
	
	h2_norm = 0;
	
	deltavect = sym(zeros(db+1,1));
	
    if(da > db)
        parfor ii=1:(db+1)
			disp(strcat('index = ',num2str(ii)));
            M_tmp   = M_hurwitz;
			idx = db - (ii-1) + 1;
            M_tmp   = simplify(setcol_an(M_tmp,poly_a,p,idx));
            deltam  = simplify(det(M_tmp));
            bm      = simplify(compute_bm(poly_b,p,ii));
			deltavect(ii) = simplify(deltam*bm);
        end
		
		h2_norm = simplify(sum(deltavect));
		
        extmp   = simplify(1/(2*a(da+1)^2*det(M_hurwitz)));
        h2_norm = simplify(h2_norm*extmp);
        extmp   = simplify(b(db + 1)*b(db)/(a(da + 1)^2));
        h2_norm = simplify(h2_norm - extmp);
		
    elseif(da == db)
		
		parfor ii=2:(db+1)
			disp(strcat('index = ',num2str(ii)));
			M_tmp   = M_hurwitz;
			idx     = db - (ii-1) + 1;
			M_tmp   = simplify(setcol_an(M_tmp,poly_a,p,idx));
			deltam  = simplify(det(M_tmp));
			bm      = simplify(compute_bm_eq(poly_a,poly_b,p,ii));
			deltavect(ii) = simplify(deltam*bm);
		end
		
		h2_norm = sum(deltavect);
		
		extmp   = simplify(1/(2*a(da+1)^2*det(M_hurwitz)));
		h2_norm = h2_norm*extmp;
		
		c1 = b(1)*(b(db)/b(1)   - a(da)/a(1)  );
		c2 = b(1)*(b(db+1)/b(1) - a(da+1)/a(1));	
		
		extmp   = simplify(c1*c2/(a(da + 1)^2));
		h2_norm = simplify(h2_norm - extmp);
		
    else
        disp('non proper transfert function cannot be treated');
	end

	h2_step = h2_norm;

	digitsOld = digits(5);
	h2_step = vpa(h2_step,4);

	[num,den] = numden(h2_step);

	h2_step = recurse_horner_new(num)/recurse_horner_new(den);
	
end
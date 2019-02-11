

%% Give stability criteria with Routh Table criteria :
%% https://en.wikipedia.org/wiki/Routh%E2%80%93Hurwitz_stability_criterion


function criteria_coefs = build_criteria_stability_routh_table(ftbf, coefs, p)

	[num_ftbf, den_ftbf] = numden(ftbf);
	
    den_deg   = symtbx_poly_degree(den_poly,p);
	den_coeff = symtbx_poly_coeffs(den_poly,p);
	
	nbcol = 0;
	if (mod(den_deg+1,2) == 1)
		nbcol = (den_deg+2)/2;
	else 
		nbcol = (den_deg+1)/2;
	end
	
    routhtable = zeros(den_deg+1,nbcol);
	
    for ii=1:2:den_deg+1
        routhtable(1,ii/2) = den_coeff(den_deg - ii);
        routhtable(2,ii/2) = den_coeff(den_deg - (ii+1));
	end
	
	criteria_coefs = {}
    coefs{end+1} = routhtable(1,1);
    coefs{end+1} = routhtable(2,1);
	
    for ii=2:den_deg+1
        for jj=1:nbcol
            if(routhtable(ii-1,jj+2) == 0 || routhtable(ii,1) == 0 )
                routhtable(ii,jj) = 0;
            else
                routhtable(ii,jj) = -1/routhtable(ii,1)*(routhtable(ii-1,1)*routhtable(ii-1,jj+1)-routhtable(ii-1,1)*routhtable(ii-2,jj+1));
                routhtable(ii,jj) = simplify(routhtable(ii,jj));
            end
        end
        criteria_coefs{end+1} = routhtable(ii,1);
    end
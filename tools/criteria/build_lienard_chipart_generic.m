function coefs = build_lienard_chipart_generic(order, name)

    for ii=1:order
        eval(['syms ' name num2str(ii) ';']);
    end

	%% Build generic hurwitz matrix
	d = order -1;
	for jj=1:d
	    even_idx = 2*(jj);
        for ii=1:d
           idx = (even_idx - ii);
            if (idx > d)
                M(ii,jj) = 0;
	    elseif (idx >= 0)
	        eval(['M(ii,jj) = ', name, num2str(idx+1),';']);
            end
        end
	end

    disp(M);
    coefs = {};
    order = order-1
    mode  = mod(order,2);

    if (mode == 0)
        ii_c = (order+1):-2:1;
        ii_m = (order)  :-2:1;
    else
        ii_c = (order+1):-2:1;
        ii_m = (order-1):-2:1;
    end

    s_c = length(ii_c);
    s_m = length(ii_m);

    for ii = 1:s_c
        idx = ii_c(ii);
        eval(['coefs{end+1} = ', name, num2str(idx),';']);
    end
	
    for ii = 1:s_m
        idx = ii_m(ii);
        minor = symtbx_minor_matrix(M,idx);
        coefs{end+1} = det(minor);
    end
    

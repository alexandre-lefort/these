%% Return a matlab transfer function from a symbolic transfert function

function tf_t = symtbx_symtf_to_tf(symtf)

    [num,den] = numden(symtf);
    
    num_c = eval(symtbx_get_sympoly_coeffs(num));
    den_c = eval(symtbx_get_sympoly_coeffs(den));
    
    if (strcmp(class(num_c),'double') ~= 1)
    	error('Remanent symbolic parameters in the expression (numerator).');
    end
    if (strcmp(class(den_c),'double') ~= 1)
    	error('Remanent symbolic parameters in the expression (denominator).');
    end

    tf_t = tf(num_c,den_c);

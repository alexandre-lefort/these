%% Return a matlab transfer function from a symbolic transfert function

function symbode = symtbx_symtf_bode(symtf, w)

    p = symvar(symtf);

    if (length(p) == 1 && p(1) ~= 's')
       error('remanent uncertain parameter');
    end

    if (length(p) > 1)
       error('remanent uncertain parameter');
    end

    tf_t  = symtbx_symtf_to_tf(symtf);
    [zz1_mag, zz1_phase, zz1_wout] = bode(tf_t, w);  
    symbode = squeeze(zz1_mag);

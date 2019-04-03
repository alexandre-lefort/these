function gabarits = build_weightings(model_sys)

	syms s;
	logw = -3:0.01:1;
	w    = 10.^(logw); %% Pulsation vector

    %% state space representation
    sym_zz0 = model_sys.G0(1,1);
    tf_zz0  = symtbx_symtf_to_tf(sym_zz0);
    [zz0_mag, ~, ~] = bode(tf_zz0,w);  
    
	%% Find cutoff frequency

	%% Build weightings

	gabarits.zz1 = s*sym(2); 
	gabarits.zz2 = sym(16/10);
	gabarits.zb = sym(16/10);

	tf_wzz1 = symtbx_symtf_to_tf(gabarits.zz1);
    tf_wzz2 = symtbx_symtf_to_tf(gabarits.zz2);

    [wzz1_mag, ~, ~] = bode(tf_wzz1, w);  
    [wzz2_mag, ~, ~] = bode(tf_wzz2, w);  

    zz0  = squeeze(zz0_mag); 
    wzz1 = squeeze(wzz1_mag);
    wzz2 = squeeze(wzz2_mag);

    figure; loglog(w, zz0,'r',w , wzz1,'b', w, wzz2, 'b');
    grid on;

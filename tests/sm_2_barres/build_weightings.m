function gabarits = build_weightings(model_sys)

	delta_w = 1.0 ; %% Pulsation margin

	syms s;
	logw = -3:0.01:1;
	w    = 10.^(logw); %% Pulsation vector

    %% state space representation
    sym_zz0 = model_sys.G0(1,1);
    tf_zz0  = symtbx_symtf_to_tf(sym_zz0);
    [zz0_mag, zz0_phase, zz0_wout] = bode(tf_zz0,w);  
    
	%% Find cutoff frequency

	cut_off  = 2;

	%% Build weightings

	gabarits.zz1 = s*sym(2); 
	gabarits.zz2 = sym(16/10);
	gabarits.zb1 = sym(16/10);
    gabarits.zb2 = sym(16/10);

	tf_wzz1 = symtbx_symtf_to_tf(gabarits.zz1);
    tf_wzz2 = symtbx_symtf_to_tf(gabarits.zz2);
    tf_wzb1 = symtbx_symtf_to_tf(gabarits.zb1);
    tf_wzb2 = symtbx_symtf_to_tf(gabarits.zb2);

    [wzz1_mag, wzz1_phase, wzz1_wout] = bode(tf_wzz1, w);  
    [wzz2_mag, wzz2_phase, wzz2_wout] = bode(tf_wzz2, w);  
    [wzb1_mag, wzb1_phase, wzb1_wout] = bode(tf_wzb1, w); 
    [wzb2_mag, wzb2_phase, wzb2_wout] = bode(tf_wzb2, w); 
    zz0  = squeeze(zz0_mag); 
    wzz1 = squeeze(wzz1_mag);
    wzz2 = squeeze(wzz2_mag);
    wzb1 = squeeze(wzb1_mag);
    wzb2 = squeeze(wzb2_mag);

    figure; loglog(w, zz0,'r',w , wzz1,'b', w, wzz2, 'b');
    grid on;

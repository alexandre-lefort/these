function criterias = build_criterias(transferts, gabarits)

    %% Build criterias separated for each variable

    criterias.zz1    = build_criteria_hinf(transferts.zz  / gabarits.zz1);
    criterias.zz2    = build_criteria_hinf(transferts.zz  / gabarits.zz2);
    criterias.zb     = build_criteria_hinf(transferts.zb  / gabarits.zb);
    criterias.zz_t   = build_criteria_hinf(transferts.zz);
    criterias.zb_t   = build_criteria_hinf(transferts.zb);
    criterias.zz1_w  = build_criteria_hinf(gabarits.zz1 );
    criterias.zz2_w  = build_criteria_hinf(gabarits.zz2 );
    criterias.zb_w   = build_criteria_hinf(gabarits.zb );
    %% Build stability criterion

    criterias.stab_coefs = build_criteria_stab_coeffs(transferts.poly_stab);
    criterias.stab_lc    = build_criteria_stab_Lienard_Chipart(transferts.poly_stab);
function criterias = build_criterias(transferts)

    %% Build stability criterion

    criterias.stab_coefs = build_criteria_stab_coeffs(transferts.poly_stab);
    criterias.stab_lc    = build_criteria_stab_Lienard_Chipart(transferts.poly_stab);

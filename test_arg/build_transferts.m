function transferts = build_transferts(model_sys, model_ctrl)

    syms s;
    
    %% Construire la fonction de transfert en immersion
    
    %% X = [dZ dTheta Z Theta]
    %% Y = [dZ dTheta Z Theta]
    %% U = [Beta1]

	K  = model_ctrl.symtf; %% controller transfert function
    G  = model_sys.G;      %% system transfert function

    %% Sensitivity
    H1 = inv(eye(4) + G*K);
    H2 = K*H1;
    
    zz  = H1(1,1);
    zb1 = H2(1,1);
    zb2 = H2(2,1);

    transferts.zz  = simplify(zz); %% Closed loop transfert function for the states
    transferts.zb1 = simplify(zb1); %% Closed loop transfert function for the commands
    transferts.zb2 = simplify(zb2); %% Closed loop transfert function for the commands
    [transferts.poly_stab,transferts.mat_stab] = symtbx_closed_loop_poly(model_sys.sys,model_ctrl.symss); %% Intern stability
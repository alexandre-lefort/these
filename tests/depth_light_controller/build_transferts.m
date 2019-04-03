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
    
    zz = H1(1,1);
    zb = H2(1,1);

    transferts.zz = simplify(zz); %% Closed loop transfert function for the states
    transferts.zb = simplify(zb); %% Closed loop transfert function for the commands
    [transferts.poly_stab,transferts.mat_stab] = symtbx_closed_loop_poly(model_sys.sys,model_ctrl.symss); %% Intern stability
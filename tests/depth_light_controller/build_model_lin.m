function sym_model = build_model_lin(Vs)
    
    model = load_model('iver2.nav');
    prec = 10;
    p = symtbx_convert_struct_to_sym(model,prec);
    p.Vs = symtbx_convert_double_to_sym(Vs,prec);
    
    ptmp = p;
    p.CzB1 = sym('CzB1');
    p.CmB1 = sym('CmB1');
    p.CzW  = sym('CzW');
    p.CmQ  = sym('CmQ');
    
    %%[Z Theta W Q beta1 0]

    Cz = sym([0 0 p.CzW/p.Vs p.CzQ*p.L/p.Vs p.CzB1 0]);
    Cm = sym([0 0 p.CmW/p.Vs p.CmQ*p.L/p.Vs p.CmB1 0]);

    Fz_navire = sym(p.Vs^2/(2*p.CRemp*p.L)*Cz);
    My_navire = sym(p.Vs^2/(2*p.CRemp)*Cm);

    Fz_couplage = sym([0 0 0       p.mu*p.Vs 0 0]);
    My_couplage = sym([0 0 0 -p.mu*p.Xg*p.Vs 0 0]);

    Fz_hydrostat = sym([0              0 0 0 0   (p.mu-1)*p.g]);
    My_hydrostat = sym([0 -p.mu*p.Zg*p.g 0 0 0 -p.mu*p.Xg*p.g]);

    Fz = sym(Fz_hydrostat + Fz_couplage + Fz_navire);
    My = sym(My_hydrostat + My_couplage + My_navire);

    Mpropre  = sym([      p.mu   -p.mu*p.Xg  ;...
                    -p.mu*p.Xg p.L^2*p.X_2]) ;
                
    Majoutee = sym([       p.Mu_3       -p.L*p.Nu_35  ;...
                     -p.L*p.Nu_35  p.L^2*p.Lambda_2]) ;
                 
    Masse = sym(Mpropre + Majoutee);

    Aux = sym(inv(Masse)*[Fz ; My]);
    
    Wpt = sym(Aux(1,:));
    Qpt = sym(Aux(2,:));
    
    thetapt = sym([0     0 0 1 0 0]);
    Zpt     = sym([0 -p.Vs 1 0 0 0]);

    F = sym([Zpt ; thetapt ; Wpt ; Qpt]);
    A = sym(F(:,1:4));
    B = sym(F(:,5));

    %% Conversion (Z Theta W Q) -> (Z Theta Pi Q)

    ZThetaPiQ = sym([1 0 0 0 ; 0 1 0 0 ; 0 1 -1/p.Vs 0 ; 0 0 0 1]); %% Transformee en X = [Theta, Pente, Q];

    A = ZThetaPiQ*A*inv(ZThetaPiQ);
    B = ZThetaPiQ*B;
    
	sys.a = simplify(A);
	sys.b = simplify(B);
	sys.c = sym(eye(4));
	sys.d = sym(zeros(4,1));
    
	%% Build central model
    
    var_po = [symvar(sys.a) symvar(sys.b)];
    p0 = [];
    for ii=1:length(var_po)
        p0.(char(var_po(ii))) = ptmp.(char(var_po(ii)));
    end
    
	sys0.a = eval(symtbx_sym_subs_from_struct(sys.a, p0));
	sys0.b = eval(symtbx_sym_subs_from_struct(sys.b, p0));
	sys0.c = eval(symtbx_sym_subs_from_struct(sys.c, p0));
	sys0.d = eval(symtbx_sym_subs_from_struct(sys.d, p0));
    
	syms s;
	
	G  = simplify(symtbx_symss_to_symtf(sys));
	G0 = simplify(symtbx_symss_to_symtf(sys0));

	%% Save in struct
	sym_model.sys   = sys  ;
	sym_model.p     = p    ;
	sym_model.G     = G    ;
    sym_model.sys0  = sys0 ;
    sym_model.p0    = p0   ;
	sym_model.G0    = G0   ;

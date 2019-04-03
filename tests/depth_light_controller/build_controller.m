%%
%% Build a depth controller using stern plane with an integrator on Z :
%%
%% X = [Z Theta Pi Q]
%% Y = [Z Theta Pi Q]
%% U = [Beta1]

function model_ctrl = build_controller(Vs,g)
     
    delay = sym('314/100');
    
    %% Controller gains 
    if (nargin == 1)
        g.kz     = sym('kz'    , 'real');
        g.ktheta = sym('ktheta', 'real');
        g.kpi    = sym('kpi'   , 'real');
        g.kq     = sym('kq'    , 'real');
        g.iz     = sym('iz'    , 'real');
    end
    
    syms s;
    adim_iz = sym('5/100');
    %ret = (2/delay)/(s +2/delay);
    
    %% Transfert matrix representation
    sym_tf = [(g.kz+adim_iz*g.iz/s) g.ktheta g.kpi g.kq];
    
    %% State space representation     
    a = zeros(1);
    b = [1 0 0 0];
    c = adim_iz*g.iz;
    d = [g.kz g.ktheta g.kpi g.kq];

    %% Add first order delay
    ret = 1/(1+s/delay);
    af = -1/delay;
    bf = -1;
    cf =  1;
    df =  0;

    sym_tf   = ret*sym_tf;
    sym_ss.a = [a 0; bf*c af];
    sym_ss.b = [b;bf*d];
    sym_ss.c = [df*c cf];
    sym_ss.d = df*d;

    %% Save struct
    model_ctrl.gains = g;
    model_ctrl.symtf = sym_tf;
    model_ctrl.symss = sym_ss;

    model_ctrl.gnames.kz     = sym('kz'    , 'real');
    model_ctrl.gnames.ktheta = sym('ktheta', 'real'); 
    model_ctrl.gnames.kpi    = sym('kpi'   , 'real'); 
    model_ctrl.gnames.kq     = sym('kq'    , 'real'); 
    model_ctrl.gnames.iz     = sym('iz'    , 'real');

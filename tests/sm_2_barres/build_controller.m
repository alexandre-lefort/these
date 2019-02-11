%%
%% Build a depth controller using stern plane with an integrator on Z :
%%
%% X = [Z Theta dZ dTheta]
%% Y = [Z Theta dZ dTheta]
%% U = [Beta1]

function model_ctrl = build_controller(gains)

    %% Controller gains 
    if (nargin == 0)
        gains.kpz1  = sym('kpz1' ,'real'); 
        gains.kdz1  = sym('kdz1' ,'real'); 
        gains.kpt1  = sym('kpt1' ,'real'); 
        gains.kdt1  = sym('kdt1' ,'real');
        gains.kiz1  = sym('kiz1' ,'real');
        gains.kpz2  = sym('kpz2' ,'real'); 
        gains.kdz2  = sym('kdz2' ,'real'); 
        gains.kpt2  = sym('kpt2' ,'real'); 
        gains.kdt2  = sym('kdt2' ,'real');
        %%gains.kiz2  = sym('kiz2' ,'real');
    end
    
    syms s;
    
    %% Transfert matrix representation
    sym_tf = [(gains.kpz1+gains.kiz1/s) gains.kpt1 gains.kdz1 gains.kdt1  ; ...
              (gains.kpz2             ) gains.kpt2 gains.kdz2 gains.kdt2 ];
    
    %% State space representation     
    sym_ss.a = zeros(1);
    sym_ss.b = [gains.kiz1 0 0 0];
    sym_ss.c = ones(2,1);
    sym_ss.d = [gains.kpz1 gains.kpt1 gains.kdz1 gains.kdt1 ; ...
                gains.kpz2 gains.kpt2 gains.kdz2 gains.kdt2];

    %% Save struct
    model_ctrl.gains = gains;
    model_ctrl.symtf = sym_tf;
    model_ctrl.symss = sym_ss;

    model_ctrl.gnames.kpz1  = sym('kpz1' ,'real'); 
    model_ctrl.gnames.kdz1  = sym('kdz1' ,'real'); 
    model_ctrl.gnames.kpt1  = sym('kpt1' ,'real'); 
    model_ctrl.gnames.kdt1  = sym('kdt1' ,'real');
    model_ctrl.gnames.kiz1  = sym('kiz1' ,'real');
    model_ctrl.gnames.kpz2  = sym('kpz2' ,'real'); 
    model_ctrl.gnames.kdz2  = sym('kdz2' ,'real'); 
    model_ctrl.gnames.kpt2  = sym('kpt2' ,'real'); 
    model_ctrl.gnames.kdt2  = sym('kdt2' ,'real');
    %%model_ctrl.gnames.kiz2  = sym('kiz2' ,'real');
	
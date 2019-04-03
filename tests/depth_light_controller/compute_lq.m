clear all
clc
close all

Vs = 1.5;
model_sm = build_model_lin(Vs);


Q = eye(4);
R = 1;
Klq = lqr(model_sm.sys0.a, model_sm.sys0.b, Q, R) ;

ctrl_gains.kz     = Klq(1);
ctrl_gains.ktheta = Klq(2);
ctrl_gains.kpi    = Klq(3);
ctrl_gains.kq     = Klq(4);
ctrl_gains.iz     = 0.0;


syms s;

logw = -2:0.0001:1;
w    = 10.^(logw); %% Pulsation vector

disp(ctrl_gains);

load('criterias.mat');

czz1 = symtbx_sym_subs_from_struct(criterias.zz1, model_sm.p0);
czz1 = simplify(symtbx_sym_subs_from_struct(czz1, ctrl_gains));
czz1_bode  = eval(max(abs(subs(czz1,'w', w))));

disp(czz1_bode);

czz2 = symtbx_sym_subs_from_struct(criterias.zz2, model_sm.p0);
czz2 = simplify(symtbx_sym_subs_from_struct(czz2, ctrl_gains));
czz2_bode  = eval(max(abs(subs(czz2,'w', w))));

disp(czz2_bode);

czb = symtbx_sym_subs_from_struct(criterias.zb, model_sm.p0);
czb = simplify(symtbx_sym_subs_from_struct(czb, ctrl_gains));
czb_bode  = eval(max(abs(subs(czb, 'w', w))));

disp(czb_bode);


for ii=1:5
        crit = criterias.stab_lc{ii};
        crit = symtbx_sym_subs_from_struct(crit, model_sm.p0);
        crit = symtbx_sym_subs_from_struct(crit, ctrl_gains);
        disp(eval(crit));
end

disp(max(max(max(czz1_bode),max(czz2_bode)),max(czb_bode)));

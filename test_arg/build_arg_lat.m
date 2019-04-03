
%% Build the function g to integrate:
poly_stab = transferts.poly_stab;
diff_poly_stab = diff(poly_stab,'s');
g = diff_poly_stab/poly_stab; 

%% Explicit real parameters in g:
p  = symvar(g);
np = length(p);
    
for ii=1:np
    scrip = strcat(char(p(ii)), ' = sym(''',char(p(ii)),''',''real'');');
    eval(scrip);
    g = subs(g,p(ii), p(ii));
end

%% Change variable to integrate on imaginary axis:

% Change s laplace variable to pulsation 1i*omega (w)
rho = sym('25');
w = sym('w','real');
g1 = subs(g, 's', 1i*w);
g2 = subs(g, 's',rho*exp(1i*w));

%% Take only real part :
g1 = real(g1);
g2 = real(g2*rho*exp(1i*w));

disp(simplify(g1))

symtbx_save_criterion(g1 , model_ctrl.gains , model_sm.p0 , 'g1.txt');
symtbx_save_criterion(g2 , model_ctrl.gains , model_sm.p0 , 'g2.txt');

%% Fix values for criterion

ctrl.gains.kz     = 1.0712  ;
ctrl.gains.ktheta = 0.9998  ;
ctrl.gains.kpi    = -10.3220;
ctrl.gains.kq     = -9.0197 ;
ctrl.gains.iz     =   0.001 ;


g1 = symtbx_sym_subs_from_struct(g1, model_sm.p0); 
g1 = symtbx_sym_subs_from_struct(g1, ctrl.gains);

g1 = simplify(g1);

[res_1,res_2] = sym_integral(g1, 0:0.1:25, 'w');
res_1 = res_1/(2*pi);
res_2 = res_2/(2*pi);

disp(res_1);
disp(res_2);

g2 = symtbx_sym_subs_from_struct(g2, model_sm.p0); 
g2 = symtbx_sym_subs_from_struct(g2, ctrl.gains);

g2 = simplify(g2);

[res_11,res_22, res] = sym_integral(g2, 0:0.1:pi/2, 'w');

res_11 = res_11/(2*pi);
res_22 = res_22/(2*pi);

disp(res_11);
disp(res_22);
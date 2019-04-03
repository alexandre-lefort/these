
%% Build the function g to integrate:
poly_stab = transferts.poly_stab;
diff_poly_stab = diff(poly_stab,'s');
g = diff_poly_stab/poly_stab; 
g_alt = poly_stab;

%% Explicit real parameters in g:
p  = symvar(g);
np = length(p);
    
for ii=1:np
    scrip = strcat(char(p(ii)), ' = sym(''',char(p(ii)),''',''real'');');
    eval(scrip);
    g = subs(g,p(ii), p(ii));
end

for ii=1:np
    scrip = strcat(char(p(ii)), ' = sym(''',char(p(ii)),''',''real'');');
    eval(scrip);
    g_alt = subs(g_alt,p(ii), p(ii));
end

%% Change variable to integrate on imaginary axis:

% Change s laplace variable to pulsation 1i*omega (w)
rho = sym('250');
w = sym('w','real');
g1 = subs(g, 's', 1i*w);
g2 = subs(g, 's',1e5*exp(1i*w));
g_alt = subs(g_alt,'s',1i*w);

diff_g_alt = diff(g_alt,'w');
g1 = diff_g_alt/g_alt;

%% Take only real part :
g1 = real(g1);
g2 = g2*1e5*exp(1i*w);

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

[res_1,res_2] = sym_integral(g1, 0:0.1:250, 'w');
res_1 = res_1/(pi);
res_2 = res_2/(pi);

disp(res_1);
disp(res_2);

g2 = symtbx_sym_subs_from_struct(g2, model_sm.p0); 
g2 = symtbx_sym_subs_from_struct(g2, ctrl.gains);

g2 = simplify(g2);

[res_11,res_22] = sym_integral(g2, -pi/2:0.1:pi/2, 'w');

res_11 = res_11/(pi);
res_22 = res_22/(pi);

disp(res_11);
disp(res_22);

g_alt = symtbx_sym_subs_from_struct(g_alt, model_sm.p0); 
g_alt = symtbx_sym_subs_from_struct(g_alt, ctrl.gains);

b = 100;
a = 0.0;
rgalt = real(g_alt);
igalt = imag(g_alt);
rdgalt = real(diff(g_alt,'w'));
idgalt = imag(diff(g_alt,'w'));

galtint1 = (rdgalt*rgalt + rdgalt*rgalt)/(rgalt*rgalt + igalt*igalt);
galtint2 = (idgalt*rgalt - rdgalt*igalt)/(rgalt*rgalt + igalt*igalt);
galtint  = simplify(galtint1 + 1i*galtint2);
disp(isequal(galtint,g1));

[resgalt1, resgalt2] = sym_integral(galtint,0:0.1:100,'w');

R1 = simplify(real(g_alt)*real(g_alt) + imag(g_alt)*imag(g_alt));
R1_b = eval(subs(R1,'w',b));
R1_a = eval(subs(R1,'w',a));

res_alt = 1/2*log(R1_b/R1_a)

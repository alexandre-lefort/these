
%% Build the function g to integrate:
poly_stab = transferts.poly_stab;
diff_poly_stab = diff(poly_stab,'s');
g = diff_poly_stab/poly_stab; 

%% Expicit real parameters in g:
p  = symvar(g);
np = length(p);
    
for ii=1:np
    scrip = strcat(char(p(ii)), ' = sym(''',char(p(ii)),''',''real'');');
    eval(scrip);
    g = subs(g,p(ii), p(ii));
end

%% Change variable to integrate on imaginary axis:

% Change s laplace variable to pulsation 1i*omega (w)
w = sym('w','real');
g1 = subs(g, 's', 1i*w);

%% Take only real part :
g1 = real(g1);

disp(simplify(g1))

symtbx_save_criterion(g1 , model_ctrl.gains , model_sm.p0 , 'g1.txt');
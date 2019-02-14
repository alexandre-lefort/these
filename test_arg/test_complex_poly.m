%%build_problem;

%% Build the function g to integrate:
g = transferts.poly_stab;

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
g = subs(g, 's', 1i*w);

real_g = simplify(real(g));
imag_g = simplify(imag(g));

diff_g = simplify(diff(g,'1i*w'));

diff_real_g = simplify(diff(real_g,'w'));
diff_imag_g = simplify(diff(imag_g,'w'));

real_diff_g = simplify(real(diff_g));
imag_diff_g = simplify(imag(diff_g));

res1 = isequal(diff_real_g,real_diff_g)
res2 = isequal(diff_imag_g,imag_diff_g)
val = simplify(1i*diff_real_g -diff_imag_g);
res3 = isequal(diff_g,val)

%% Unit Test symtbx_sym_subs_from_struct

addpath('..');

syms a b c;
p.a  = 1;
m.a  = c;
expr = a+b+c;

res1 = symtbx_sym_subs_from_struct(expr, p);
res2 = symtbx_sym_subs_from_struct(expr, m);

res1_sol = 1+b+c;
res2_sol = b+2*c;

assert(isequal(res1,res1_sol) == true);
assert(isequal(res2,res2_sol) == true);

addpath('../symtbx');

omega = 0.51;
z     = -1.1;

s = tf('s');

eq = omega^2/(s^3 + s^2 + 2*z*omega*s + omega^2);

step(eq);

syms p;

eq2 =  omega^2/(p^3 + p^2 + 2*z*omega*p + omega^2);

pole(eq)

den = (p^3 + p^2 + 2*z*omega*p + omega^2);


criterias_coefs = build_criteria_stab_hurwitz_det(den, p);

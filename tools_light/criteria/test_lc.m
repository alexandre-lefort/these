syms s;
syms a0 a1 a2 a3 a4 a5 a6 a7;

poly = a0*s^7 + a1*s^6 + a2*s^5 + a3*s^4 + a4*s^3 + a5*s^2 + a6*s + a7;

c = build_criteria_stab_Lienard_Chipart(poly,s);

for ii=1:length(c)
c{ii}
horner(c{ii})
end
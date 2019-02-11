clear all
clc 
close all 

addpath('../sources');

syms p;

syms a0 a1 a2 a3 a4 a5;
den = a0*p^5 + a1*p^4 + a2*p^3 + a3*p^2 + a4*p + a5;

criteria_coefs = build_criteria_stab_Lienard_Chipart(den, p);

for ii=1:length(criteria_coefs)
    disp(simplify(criteria_coefs{ii}));
end
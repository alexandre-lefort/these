clear all
clc
close all


syms x1 x2 x3 a1 a2 a3
p = 3*a1*x1^4*(x2^2 - 8*x1*x2^2*x3 + 3*x1^5)*x3 + 6*x1^6 + a2*a3*x3^2*x2 + a1*x1*x2 + a2*4*x1^2*x2 + 6*x3 + a3 + 4;
%% p = x1^6 + x2^2 + x1^2*x2*3 -5*x1 -a3 + 5;
%% p = 6*x1^6 + 3*a1*x1^4*x2^2 - 8*x1*x2^2*x3 + a2*a3*x3^2*x2 + a1*x1*x2 + a2*4*x1^2*x2 + 3*x1^5*x3 + 6*x3 + a3 + 4;
p

x1_t = -2:2:2;
x2_t = -2:2:2;
x3_t = -2:2:2;
a1_t = -2:2:2;
a2_t = -2:2:2;
a3_t = -2:2:2;

p_new = recurse_horner_new(p);

p
p_new

p_exp = expand(p_new);

p_exp

for ii=1:3
ii
for jj=1:3
for kk=1:3
for ll=1:3
for mm=1:3
for nn=1:3

    val_p     = subs(p    ,{'x1' 'x2' 'x3' 'a1' 'a2' 'a3'},{x1_t(ii) x2_t(jj) x3_t(kk) a1_t(ll) a2_t(mm) a3_t(nn)});
    val_p_new = subs(p_new,{'x1' 'x2' 'x3' 'a1' 'a2' 'a3'},{x1_t(ii) x2_t(jj) x3_t(kk) a1_t(ll) a2_t(mm) a3_t(nn)});

    %% val_p
    %% val_p_new

    if (val_p ~= val_p_new)
        error('failed');
    end

end
end
end
end
end
end


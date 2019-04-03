clear all
clc 
close all 

addpath('../sources');

s = tf('s');
syms p;

for a = -0.5:0.2:0.5

            a
            

    for     b = -1.5:0.2:1.5
        for c = -3.0:0.4:3.0


            eq = c/(s^3 + s^2 + 2*b*c*s + a^2);

            stable = 1;

            if max(real(pole(eq))) >=0
                stable = 0;
            end

            den = (p^3 + p^2 + 2*b*c*p + a^2);

            criteria_coefs = build_criteria_stab_Lienard_Chipart(den, p);

            valid = 1;
            for ii=1:length(criteria_coefs)
                if (criteria_coefs{ii}<=0)
                    valid = 0;
                end
            end

            valid
            if valid ~= stable
                eq
                error('bug');
            end

        end
    end
end


%% Return the hinf criteria of a problem as its modulus
%%

function criteria_hinf = build_criteria_hinf(mod)

    %% Define all parameters as real values

    p  = symvar(mod);
    np = length(p);

    for ii=1:np
        scrip = strcat(char(p(ii)), ' = sym(''',char(p(ii)),''',''real'');');
        eval(scrip);
        mod = subs(mod,p(ii), p(ii));
    end

    %% Change s laplace variable to pulsation 1i*omega (w)

    w = sym('w','real');
    mod = subs(mod, 's', 1i*w);
	
    %% Square module of the transfert function as hinf criteria

   [num_mod,den_mod] = numden(mod);

    numD = real(num_mod)^2 + imag(num_mod)^2;
    denD = real(den_mod)^2 + imag(den_mod)^2;

    criteria_hinf = numD/denD;

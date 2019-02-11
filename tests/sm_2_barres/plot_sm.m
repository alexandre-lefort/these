function plot_sm(model_sm, p0vect, size_f)

    syms s;

    logw = -3:0.01:1;
    w    = 10.^(logw); %% Pulsation vector

    nCzW = length(p0vect.CzW);
    nCmQ = length(p0vect.CmQ);
    
    %% Transfert to Z
    sym_zz = model_sm.G(1,1);

    %% Bode Z

    fg =figure('units','normalized','outerposition',[0 0 1 1]);hold on;

    x = xlabel('Pulsation (rad/s)');
    set(x,'Interpreter','latex','FontSize',size_f);
    y = ylabel('Magnitude (dB)');
    set(y,'Interpreter','latex','FontSize',size_f);

    %% Plot closed-loop transfer functions

    for ii=1:nCzW
    for ll=1:nCmQ

        model_sm.p0.CzW = p0vect.CzW(ii) ;
        model_sm.p0.CmQ = p0vect.CmQ(ll) ;

        sym_zz0 = symtbx_sym_subs_from_struct(sym_zz, model_sm.p0); 
        zz0     = symtbx_symtf_bode(sym_zz0); 

        hold on;
        set(gca, 'YScale', 'log')
        set(gca, 'XScale', 'log')
        loglog(w, zz0,'b'); grid on;

    end
    end

    h = legend('$T_{b \rightarrow z}(p^*)$');
    set(h,'Interpreter','latex','FontSize',size_f);

    saveas(fg,'sm_transfert.png');
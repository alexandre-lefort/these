function plot_bodes(criterias, ctrl, model_sm, p0vect,size_f)

    syms s;

    logw = -3:0.01:1;
    w    = 10.^(logw); %% Pulsation vector

    nCzW = length(p0vect.CzW);
    nCmQ = length(p0vect.CmQ);
    
    %% Bode Z

    fg1 = figure('units','normalized','outerposition',[0 0 1 1]);hold on;

    x = xlabel('Pulsation (rad/s)'); set(x,'Interpreter', 'latex', 'FontSize', size_f);
    y = ylabel('Magnitude (dB)'   ); set(y,'Interpreter', 'latex', 'FontSize', size_f);
    
    %% Plot weighting function

    zz1_w   = symtbx_symtf_permut_s_w(criterias.zz1_w);
    zz2_w   = symtbx_symtf_permut_s_w(criterias.zz2_w);

    zz1 = symtbx_symtf_bode(zz1_w, w);
    zz2 = symtbx_symtf_bode(zz2_w, w);

    zz_w  =  min(zz1, zz2);

    set(gca, 'YScale', 'log');
    set(gca, 'XScale', 'log');
    loglog(w, zz_w, 'k'); grid on;    

    %% Plot closed-loop transfer functions

    for ii=1:nCzW
    for ll=1:nCmQ

        model_sm.p0.CzW = p0vect.CzW(ii) ;
        model_sm.p0.CmQ = p0vect.CmQ(ll) ;
        
        czz = symtbx_symtf_permut_s_w(criterias.zz_t);
        czz = symtbx_sym_subs_from_struct(czz, model_sm.p0); 
        czz = symtbx_sym_subs_from_struct(czz, ctrl.gains); 

        czz_bode  = symtbx_symtf_bode(czz, w);

        set(gca, 'YScale', 'log');
        set(gca, 'XScale', 'log');
        hold on; loglog(w, czz_bode, 'b');

    end
    end

    %%xlim([0.003 1]);
    %%ylim([1e-5,10]);

    h = legend('$W_{e_z \rightarrow z}^{-1}$','$T_{e_z \rightarrow z}(p^*,k^*)$','Location','northwest');
    set(h,'Interpreter','latex','FontSize',size_f);

    saveas(fg1,'bode_cl_30.png');

    %% Bode dS

    fg2 = figure('units','normalized','outerposition',[0 0 1 1]);hold on;
    x = xlabel('Pulsation (rad/s)');
    set(x,'Interpreter','latex','FontSize',size_f);
    y = ylabel('Magnitude (dB)');
    set(y,'Interpreter','latex','FontSize',size_f);
    grid on;

    %% Plot weighting function
    zb_w      = symtbx_symtf_permut_s_w(criterias.zb_w);
    zb_w_bode = symtbx_symtf_bode(zb_w, w);

    set(gca, 'YScale', 'log');
    set(gca, 'XScale', 'log');
    loglog(w, zb_w_bode, 'k'); grid on;    

    %% Plot closed-loop transfer functions

    for ii=1:nCzW
    for ll=1:nCmQ

        model_sm.p0.CzW = p0vect.CzW(ii) ;
        model_sm.p0.CmQ = p0vect.CmQ(ll) ;
        
        czb = symtbx_symtf_permut_s_w(criterias.zb_t);
        czb = symtbx_sym_subs_from_struct(czb, model_sm.p0); 
        czb = symtbx_sym_subs_from_struct(czb, ctrl.gains); 

        czb_bode  = symtbx_symtf_bode(czb, w);

        set(gca, 'YScale', 'log');
        set(gca, 'XScale', 'log');
        hold on; loglog(w, czb_bode, 'b');

    end
    end
    
    h = legend('$W_{e_z \rightarrow \delta_S}^{-1}$','$T_{e_z \rightarrow \delta_S}(p^*,k^*)$','Location','southeast');
    set(h,'Interpreter','latex','FontSize',size_f);

    saveas(fg2,'bode_cl_30_ds.png');
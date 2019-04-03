fg = figure('units','normalized','outerposition',[0 0 1 1]); hold on;

nCzW = length(p0vect.CzW);
nCmQ = length(p0vect.CmQ);

%% Initial State
init_sm      = [0;0;0;0];
initk        = 0;
ZCo_init     = 0;
ThetaCo_init = 0;
ZCo_step     = 1;
ThetaCo_step = 0;
DT           = 0.1;
T            = [0 100];

%% Controller
Ak = model_ctrl.symss.a;
Bk = model_ctrl.symss.b;
Ck = model_ctrl.symss.c;
Dk = model_ctrl.symss.d;

%% Start loop
for ii=1:nCzW
    for ll=1:nCmQ
        model_sm.p0.CzW = p0vect.CzW(ii) ;
        model_sm.p0.CmQ = p0vect.CmQ(ll) ;

        A = eval(symtbx_sym_subs_from_struct(model_sm.sys.a, model_sm.p0));
        B = eval(symtbx_sym_subs_from_struct(model_sm.sys.b, model_sm.p0));
        C = eval(model_sm.sys.c);
        D = eval(model_sm.sys.d);

        sim('sm',T);
    
        tt    = mesure.Time;
        Z     = mesure.Data(:,1);
        Theta = mesure.Data(:,2);
        Beta1 = commande.Data(:,1);

        grid on; hold on;
        
        subplot(3,1,1); hold on;
        plot(tt, Z, 'b');
        y = ylabel('Depth (m)');
        set(y,'Interpreter','latex','FontSize',size_f);
        grid on;
        
        subplot(3,1,2); hold on;
        plot(tt, Theta*180/pi, 'b');
        y = ylabel('Pitch ($^{\circ}$)');
        set(y,'Interpreter','latex','FontSize',size_f);
        grid on;
        
        subplot(3,1,3); hold on;
        plot(tt, Beta1*180/pi, 'b');
        y = ylabel('Plane ($^{\circ}$)');
        set(y,'Interpreter','latex','FontSize',size_f);
        grid on;

        x = xlabel('Time ($^{\circ}$)');
        set(x,'Interpreter','latex','FontSize',size_f);

    end
end

saveas(fg,'sm_step.png');
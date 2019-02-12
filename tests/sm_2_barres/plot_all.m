clear all
clc
close all


load('criterias.mat');

%% Selected controller

model_sm = build_model_lin();

ctrl_gains.kpz1= 1.0712    ;    
ctrl_gains.kdz1= 0.9998    ; 
ctrl_gains.kpt1= -10.3220  ;   
ctrl_gains.kdt1= -9.0197   ;  
ctrl_gains.kpz2= -0.3585   ;  
ctrl_gains.kdz2= -0.0221   ;  
ctrl_gains.kpt2= -1.8123   ;  
ctrl_gains.kdt2= -0.8398   ;  
ctrl_gains.kiz1=   0.001   ;  


model_ctrl = build_controller(ctrl_gains);

logw = -2:0.1:1.0;
w = 10.^(logw);
nw = length(w);

%% Uncertain parameters

pf = fields(model_sm.p0);

for ii=1:length(pf)
	pval   = model_sm.p0.(pf{ii});
    p0vect.(pf{ii}) = (-0.1*abs(pval):0.03*abs(pval):0.1*abs(pval)) + pval ;
end

%% Plot step response

size_f = 24;

plot_step;

%% Plot bodes

%% plot_bodes(criterias, model_ctrl, model_sm, p0vect,size_f);

%%plot_sm(model_sm, p0vect, size_f);

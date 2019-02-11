clear all
clc
close all


load('criterias.mat');

%% Selected controller

model_sm = build_model_lin();


ctrl_gains.kpz =  1.95703  ;
ctrl_gains.kdz = -3.78041  ;
ctrl_gains.kpt = -6.84096  ;
ctrl_gains.kdt = -0.487614 ;
ctrl_gains.kiz =  0.264851 ;

model_ctrl = build_controller(ctrl_gains);

logw = -2:0.01:1.0;
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

%% plot_step;

%% Plot bodes

plot_bodes(criterias, model_ctrl, model_sm, p0vect,size_f);

%% plot_sm(model_sm, p0vect, size_f);

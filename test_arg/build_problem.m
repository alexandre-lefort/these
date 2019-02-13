clear all
clc
close all

%% Build linear model
	disp('--> build linear model');    
	model_sm = build_model_lin();

%% Build controller model
	disp('--> build controller model');    
	model_ctrl = build_controller();

%% Build transfert functions
	disp('--> build transferts'); 
	transferts = build_transferts(model_sm, model_ctrl);

%% Build criterias
	disp('--> build criterias'); 	
	criterias = build_criterias(transferts);

	save('criterias.mat', 'criterias');

symtbx_save_criterion(criterias.stab_coefs , model_ctrl.gains , model_sm.p0 , 'Tstab_coefs.txt');
symtbx_save_criterion(criterias.stab_lc    , model_ctrl.gains , model_sm.p0 , 'Tstab_lc.txt'   );

symtbx_save_criterion(criterias.stab_coefs, model_ctrl.gains , model_sm.p0 , 'tstab.h');
save_lienard_chipart(length(criterias.stab_coefs), 'stab' , 'lc.h', sym(0));

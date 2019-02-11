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
	disp('--> build weightings'); 
	gabarits = build_weightings(model_sm);

%% Build criterias
	disp('--> build criterias'); 	
	criterias = build_criterias(transferts, gabarits, model_sm.p0, model_ctrl.gains);

	save('criterias.mat', 'criterias');

%% Save in files

symtbx_save_criterion(criterias.zz1        , model_ctrl.gains , model_sm.p0 , 'functions/Tzz1.txt'       );
symtbx_save_criterion(criterias.zz2        , model_ctrl.gains , model_sm.p0 , 'functions/Tzz2.txt'       );
symtbx_save_criterion(criterias.zb1        , model_ctrl.gains , model_sm.p0 , 'functions/Tzb1.txt'       );
symtbx_save_criterion(criterias.zb2        , model_ctrl.gains , model_sm.p0 , 'functions/Tzb2.txt'       );
symtbx_save_criterion(criterias.stab_coefs , model_ctrl.gains , model_sm.p0 , 'functions/Tstab_coefs.txt');
symtbx_save_criterion(criterias.stab_lc    , model_ctrl.gains , model_sm.p0 , 'functions/Tstab_lc.txt'   );

symtbx_save_criterion(criterias.stab_coefs, model_ctrl.gains , model_sm.p0 , 'ibex/tstab.h');
save_lienard_chipart(length(criterias.stab_coefs), 'stab' , sym(0));
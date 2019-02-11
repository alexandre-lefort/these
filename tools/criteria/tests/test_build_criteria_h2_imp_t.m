function [res,h2_imp_res] = test_build_criteria_h2_imp

	close all
	%% clear all
	clc

	addpath('../symtbx')
	syms a0 a1 a2 a3 b0 b1 b2;
	syms s;
	
	a0_t = 1;
	a1_t = 2;
	a2_t = 3;
	a3_t = 4;
	b0_t = 0;
    b1_t = 1;
    b2_t = 2;

	system_sym = (b0*s^2 + b1*s + b2)/(a0*s^3 + a1*s^2 + a2*s + a3);
	system_tes = subs(system_sym,{a0 a1 a2 a3 b0 b1 b2},{a0_t a1_t a2_t a3_t b0_t b1_t b2_t});


	system_num = tf([b0_t b1_t b2_t],[a0_t a1_t a2_t a3_t]);
	
	system_num
	system_tes
	
	dt = 0.00001;
	Tvect = 0:dt:100;
	[Y ,T ] = impulse(system_num,Tvect);

	plot(T,Y);

	h2_imp = build_criteria_h2_imp(system_sym, s);
	
	h2_imp
	
	res = sum((Y).*(Y)).*dt;

	h2_imp_res = eval(subs(h2_imp,{a0 a1 a2 a3 b0 b1 b2},{a0_t a1_t a2_t a3_t b0_t b1_t b2_t}));

	res
	
	h2_imp_res

	bool = 1;
	
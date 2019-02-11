function bool = test_build_criteria_h2_step

	close all
	clear all
	clc

	syms pa pb pc;
	syms s;
	
	pa_t = 1;
	pb_t = 5;
	pc_t = 8;

	G = (s+1)/(s^2 + 8*s + 1);
	K = 5 + 0.5/s;
	
	sysbo = K*G/(1+K*G);
	system_sym = sysbo;
	
	system_num = symtbx_sym2tf(sysbo,s);
	
	system_num
	system_sym
	
	dt = 0.00001;
	Tvect = 0:dt:100;
	
	[Y,T] = step(system_num,Tvect);
	
	plot(T,Y);
	
	h2_step = build_criteria_h2_step(system_sym, s);
	
	h2_step
	
	res = sum((Y-1).*(Y-1).*dt);

	h2_step_res = eval(subs(h2_step,{pa pb pc},{pa_t pb_t pc_t}));
	
	res
	
	h2_step_res
	
	bool = 1;
	
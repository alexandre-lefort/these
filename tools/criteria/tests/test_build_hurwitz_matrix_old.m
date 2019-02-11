function bool = test_build_hurwitz_matrix_old
	
	syms  p;

	%% Test for n = 5
	
	a = sym('a',[6 1]);
	
	poly = a(1)*p^5 + a(2)*p^4 + a(3)*p^3 + a(4)*p^2 + a(5)*p + a(6);
	
	M = build_hurwitz_matrix_old(poly, p);
	
	M_test = [a(2) a(4) a(6)    0    0 ; ...
			  a(1) a(3) a(5)    0    0 ; ...
				0  a(2) a(4) a(6)    0 ; ...
				0  a(1) a(3) a(5)    0 ; ...
				0     0 a(2) a(4) a(6) ];
	
	bool = isequaln(M,M_test);
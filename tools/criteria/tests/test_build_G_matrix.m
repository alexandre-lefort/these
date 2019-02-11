function bool = test_build_G_matrix

	a = sym('a',4);
	g = sym('g',[6 1]);
	
	syms p;
	
	poly = g(1)*p^5 + g(2)*p^4 + g(3)*p^3 + g(4)*p^2 + g(5)*p + g(6);
	
	M = sym([1 2 3  ; ...
			 5 6 1  ; ...
			 9 9 0]); 
	
	G_res = build_G_matrix(M, poly, p);
	
	G_test = [g(1) g(2) g(3) ; ...
				5    6    1 ; ...
				9    9    0];
				
	bool = isequaln(G_res, G_test);
function bool = test_compute_bm
	
	syms p;
	
	b = sym('b',[6 1]);
	
	poly = b(1)*p^5 + b(2)*p^4 + b(3)*p^3 + b(4)*p^2 + b(5)*p + b(6);


	B1_test = sym(zeros(6,1));
	
	B1_test(6) = b(6)^2;	
	B1_test(5) = b(5)^2 - 2*b(6)*b(4);		
	B1_test(4) = b(4)^2 - 2*b(5)*b(3) + 2*b(6)*b(2);		
	B1_test(3) = b(3)^2 - 2*b(4)*b(2) + 2*b(5)*b(1);		
	B1_test(2) = b(2)^2 - 2*b(3)*b(1);		
	B1_test(1) = b(1)^2;
	
	for ii = 1:6
		B1_res(ii,1) = compute_bm(poly, p, ii);
	end
	
	bool = isequaln(B1_res,B1_test);
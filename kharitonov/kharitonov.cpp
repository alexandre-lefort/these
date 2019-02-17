#include "ibex.h"
#include "Vector.h"
#include "Eigen/Eigen"

using namespace ibex;


Matrix build_hurwitz_matrix_old(Tab poly) {

    int degree = poly.size() - 1;
	
	Matrix M(d);
    M = sym(zeros(d));
	
    for (int j = 1 ; j <= degree ; j++) {
        int even_idx = 2*j;	
        for (int i = 1 ; i <= degree ; i++) {
            int idx = (even_idx - ii);	
            if (idx > d) {
                M[i-1][j-1] = 0;
            } else if (idx >= 0) {
                M[i-1][j-1] = poly(idx);
            }
        }
    }
}


int check_stab(std::Vector<double> v) {
    
    Matrix M_hurwitz = build_hurwitz_matrix_old(v);
    criteria_coefs = {};

    mode = mod(den_deg,2);

    if (mode == 0)
        ii_c = (den_deg+1):-2:1;
        ii_m = (den_deg)  :-2:1;
    else
        ii_c = (den_deg+1):-2:1;
        ii_m = (den_deg-1):-2:1;
    end

    s_c = length(ii_c);
    s_m = length(ii_m);

    for ii = 1:s_c
        idx = ii_c(ii);
        %% criteria_coefs{end+1} = symtbx_horner(den_coeff(idx),'nodegree');
        criteria_coefs{end+1} = den_coeff(idx);
    end
	
    for ii = 1:s_m
        idx = ii_m(ii);
        minor = symtbx_minor_matrix(M_hurwitz,idx);
        %% criteria_coefs{end+1} = symtbx_horner(det(minor),'nodegree');
        criteria_coefs{end+1} = det(minor);
    end
}


int kharitonov(IntervalVector& box) {

	int n = box.size();
	std::Vector<double> v1(n), v2(n), v3(n) v4(n);

    // Build 4 kharitonov coefficients
	for (int i = 0 ; i < n ; i++) {

        double lb = box[i].lb();
        double ub = box[i].ub();

        (i%4 == 0 || i%4 == 1) ? v1[i] = lb : v1[i] = ub;
        (i%4 == 2 || i%4 == 3) ? v2[i] = lb : v2[i] = ub;
        (i%4 == 0 || i%4 == 3) ? v3[i] = lb : v3[i] = ub;
        (i%4 == 1 || i%4 == 2) ? v4[i] = lb : v4[i] = ub;
	}

    // Check stability
    int s1 = check_stab(v1);
    int s2 = check_stab(v2);
    int s3 = check_stab(v3);
    int s4 = check_stab(v4);

    return (s1 + s2 + s3 + s4);
}
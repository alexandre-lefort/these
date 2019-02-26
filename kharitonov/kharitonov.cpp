#include "kharitonov.h"


using namespace ibex;


Eigen::MatrixXd build_hurwitz_matrix_old(Eigen::VectorXd& poly) {

    int d = poly.size() - 1;
	
	Eigen::MatrixXd M = Eigen::MatrixXd::Zero(d,d);
	
    for (int j = 1 ; j <= d ; j++) {
        int even_idx = 2*j;	
        for (int i = 1 ; i <= d ; i++) {
            int idx = (even_idx - ii);	
            if (idx > d) {
                M[i-1][j-1] = 0;
            } else if (idx >= 0) {
                M[i-1][j-1] = poly(idx);
            }
        }
    }

    return M;
}


int check_stability(Eigen::VectorXd& v) {
    
    Eigen::MatrixXd M_hurwitz = build_hurwitz_matrix_old(v);
	int n = v.size();
    bool stable = true;
    for (int i = 1 ; i < n ; i ++) {
        Eigen::MatrixXd minor = M_hurwitz.block(0,0,i,i);
        stable &= (minor.determinant() > 0);
    }
    return stable;
}


int kharitonov(Ibex::IntervalVector& box) {

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
    int s1 = check_stability(v1);
    int s2 = check_stability(v2);
    int s3 = check_stability(v3);
    int s4 = check_stability(v4);

    return (s1 + s2 + s3 + s4);
}
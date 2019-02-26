#include "dabbene.h"


// dabbene : return true if there is a stable polynomial in the box. return false if not
bool dabbene(ibex::IntervalVector& poly, int max_iter) {

    bool found = false;
    int n = poly.size();
    if (n < 2) {
        return true;
    } else {
        ibex::IntervalVector odd_coeffs  = get_odd_coeffs (poly);
        ibex::IntervalVector even_coeffs = get_even_coeffs(poly);
    
        int no = odd_coeffs.size()  - 1;
        int ne = even_coeffs.size() - 1;
    
        Eigen::MatrixXd ko, ko_roots, ve;
        int iter = 0;


        while (iter < max_iter && !found) {
            ko       = get_ko(odd_coeffs, max_iter);
            ko_roots = compute_roots_odd_polynomial(ko);
            ve       = compute_ve(ko_roots, ne);    
            found    = does_ke_exist(even_coeffs, ve);
            iter++;
        }
    }
    return found;
}


ibex::IntervalVector get_odd_coeffs(ibex::IntervalVector& poly) {

    int n = poly.size()-1; // n is >= 2
    int no = ((n%2 == 0) ? n/2 -1 : (n-1)/2);
    ibex::IntervalVector odd_coeffs(no+1);
    for (int i = 0 ; i <= no ; i++) {
        odd_coeffs[i] = poly[2*i+1];
    }
    return odd_coeffs;
}


ibex::IntervalVector get_even_coeffs(ibex::IntervalVector& poly) {

    int n = poly.size()-1; // n is >= 2
    int ne = ((n%2 == 0) ? n/2 : (n-1)/2);
    ibex::IntervalVector even_coeffs(ne+1);
    for (int i = 0 ; i <= ne ; i++) {
        even_coeffs[i] = poly[2*i];
    }
    return even_coeffs;
}


Eigen::MatrixXd get_ko(ibex::IntervalVector& odd_coeffs, int max_iter) {

    int no = odd_coeffs.size()-1;
    int iter = 0;
    bool found = false;
    ibex::IntervalVector aux0(1), aux1(1);
    Eigen::MatrixXd ko(no+1,1);
    
    while (iter < max_iter && !found) {
    	
        // Choose k1 and k3    
        aux0[0] = odd_coeffs[0];
        ko(0,0) = aux0.random()[0];
        if (no == 0) return ko;

        aux1[0] = odd_coeffs[1];
        ko(1,0) = aux1.random()[0];
        if (no == 1) return ko;

        //std::cout << iter << " " << ko(0,0) << " " << ko(1,0) << std::endl;
        
        int i = 2;
        while (i < no + 1 && iter < max_iter) {
            if(is_interval_empty(odd_coeffs[i], ko(i-1,0), ko(i-2,0), (double) i, (double) no)) {
            	//std::cout << "true " << i << std::endl;
                iter++;
                break;
            } else {
            	//std::cout << "false " << i << std::endl;
                aux1[0] = odd_coeffs[i];
                ko(i,0) = aux1.random()[0];
                i++;
            }
        }
        //std::cout << "tada " << i << " " << no+1 << std::endl;
        if (i == (no+1)) {
            found = true;
        }
    }

    return ko;
}


bool is_interval_empty(ibex::Interval k2ip1, double k2im1, double k2im3, double i, double no) {
    double c = (i-1)/i*((no-i+1)/(no-i+2));
    double lb = k2ip1.lb();
    double ub = std::min(k2ip1.ub(), c*(k2im1*k2im1)/k2im3);
    //std::cout << "lb = " << lb << " : ub = " << ub <<  " : c = " << c << std::endl;
    return (lb > ub);
}


Eigen::MatrixXd compute_roots_odd_polynomial(Eigen::MatrixXd& odd_coeffs) {

    // Build companion matrix
    //std::cout << "res 1" << std::endl;
    int no = odd_coeffs.rows()-1;
    Eigen::MatrixXd companion = Eigen::MatrixXd::Zero(no,no);
    for (int i = 0; i < no-1 ; i++) {
      companion(i+1,i) = 1;
    }
    companion.block(0,no-1,no,1) = -odd_coeffs.block(1,0,no,1)/odd_coeffs(0,0);
    // return eigenvalues as roots of the polynomial
    Eigen::VectorXcd eigenvals = companion.eigenvalues();
    Eigen::VectorXd  real_eigenvals = eigenvals.real();
    Eigen::MatrixXd  res(no,1);
    res.col(0) = real_eigenvals.col(0);
    return res;
}
   
   
Eigen::MatrixXd compute_ve(Eigen::MatrixXd& odd_roots, int ne) {
    
    int no = odd_roots.size();
    Eigen::MatrixXd ve(no+1, ne+1);
    double omega = 0;
    for (int i = 0 ; i < no+1 ; i++) {
    	if (i > 0) omega = odd_roots(i-1,0);
        Eigen::MatrixXd line = -pow(-1,i)*compute_ue(omega, ne);
        ve.block(i,0,1,ne+1) = line.block(0,0,1,ne+1);
    }
    return ve;
}


Eigen::MatrixXd compute_ue(double omega, int ne) {

    Eigen::MatrixXd ue(1,ne+1);
    for (int i = 0 ; i < ne+1 ; i++) {
        ue(0,i) = pow(-1,i)*pow(omega,2*i);
    }
    return ue;
}


bool does_ke_exist(ibex::IntervalVector& even_coeffs, Eigen::MatrixXd ve) {
    int n = ve.cols(); // no+1
    int m = ve.rows(); // ne+1
    Eigen::MatrixXd A(n+2*m, m);
    Eigen::VectorXd B(n+2*m   );
    Eigen::VectorXd c(1);
    Eigen::VectorXd f(m);
    for (int i = 0 ; i < m ; i ++) {
        B(i    ) =  0                  ;
        B(i+m  ) =  even_coeffs[i].ub();
        B(i+m*2) = -even_coeffs[i].lb();
        f(i    ) = 1;
    }
    //std::cout << "test 1 " << std::endl;
    A.block(0   ,0,n,m) = ve;
    A.block(n   ,0,m,m) = Eigen::MatrixXd::Identity(m,m);
    A.block(n+m ,0,m,m) = Eigen::MatrixXd::Identity(m,m);
    int k = n+2*m;

    return igl::linprog( c, A, B, k, f);
}
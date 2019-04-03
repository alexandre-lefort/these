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
    
        Eigen::MatrixXd roots, ve, po;
        int iter = 0;


        while (iter < max_iter && !found) {
            Eigen::MatrixXd ko(no+1,1);
            if(get_ko(ko, odd_coeffs, max_iter)) {
                //std::cout << "found ko : " << std::endl; std::cout << ko << std::endl;
                po    = compute_po(ko);
                roots = compute_roots_odd_polynomial(po);
                ve       = compute_ve(roots, ne);
                //std::cout << "ve : " << std::endl; std::cout << ve << std::endl;
                found    = does_ke_exist(even_coeffs, ve);
                //std::cout << "found ke ? : " << std::endl; std::cout << found << std::endl;
            }
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


 bool get_ko(Eigen::MatrixXd& res, ibex::IntervalVector& odd_coeffs, int max_iter) {

    int no = odd_coeffs.size()-1;
    int iter = 0;
    bool found = false;
    ibex::IntervalVector aux0(1), aux1(1);
    Eigen::MatrixXd ko(no+1,1);
    
    while (iter < max_iter && !found) {
    	
        // Choose k1 and k3    
        aux0[0] = odd_coeffs[0];
        ko(0,0) = aux0.random()[0];
        if (no == 0) {res = ko; return true;}

        aux1[0] = odd_coeffs[1];
        ko(1,0) = aux1.random()[0];
        if (no == 1) {res = ko; return true;}

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

    res = ko;
    //std::cout << found << std::endl;
    return (found);
}


bool is_interval_empty(ibex::Interval k2ip1, double k2im1, double k2im3, double i, double no) {
    double c = (i-1)/i*((no-i+1)/(no-i+2));
    double lb = k2ip1.lb();
    double ub = std::min(k2ip1.ub(), c*(k2im1*k2im1)/k2im3);
    //std::cout << "lb = " << lb << " : ub = " << ub <<  " : c = " << c << std::endl;
    return (lb > ub);
}


Eigen::MatrixXd compute_po(Eigen::MatrixXd& ko) {
    int no = ko.rows() - 1;
    Eigen::MatrixXd u  = compute_u(1, no);
    Eigen::MatrixXd po = Eigen::MatrixXd::Zero(2*no+1, 1);
    for (int i = 0 ; i < no+1 ; i++) {
        po(2*i, 0) = u(0,i)*ko(i, 0);
    }
    return (po);
}


Eigen::MatrixXd compute_roots_odd_polynomial(Eigen::MatrixXd& p) {
    // Build companion matrix
    int degree = p.rows()-1;
    Eigen::MatrixXd companion = Eigen::MatrixXd::Zero(degree, degree);
    for (int i = 0 ; i < degree-1 ; i++) {
        companion(i+1, i) = 1;
    }

    for (int i = 0 ; i < degree ; i++) {
        companion(i, degree-1) = -p(i, 0)/p(degree, 0);
    }

    Eigen::VectorXd eigenvals = companion.eigenvalues().real();
    Eigen::VectorXd eigenvals_sorted;
    Eigen::VectorXi aux;
    igl::sort(eigenvals, 1, true, eigenvals_sorted, aux);

    assert(eigenvals_sorted.size()%2 == 0);
    int neigen = eigenvals_sorted.size()/2;

    //std::cout << eigenvals_sorted << std::endl << std::endl;
    Eigen::VectorXd res(neigen);
    for (int i = 0 ; i < neigen ; i ++) {
        res(i) = eigenvals(2*i);
    }

    //std::cout << res << std::endl;
    //std::cout << "polynomial : " << std::endl; std::cout << p         << std::endl;
    //std::cout << "companion  : " << std::endl; std::cout << companion << std::endl;
    //std::cout << "roots      : " << std::endl; std::cout << eigenvals << std::endl;
    return res;
}
   
   
Eigen::MatrixXd compute_ve(Eigen::MatrixXd& odd_roots, int ne) {
    
    int no = odd_roots.size();
    Eigen::MatrixXd ve(no+1, ne+1);
    double omega = 0;
    for (int i = 0 ; i < no+1 ; i++) {
    	if (i > 0) omega = odd_roots(i-1,0);
        Eigen::MatrixXd line = -pow(-1,i)*compute_u(omega, ne);
        ve.block(i,0,1,ne+1) = line.block(0,0,1,ne+1);
    }
    return ve;
}


Eigen::MatrixXd compute_u(double omega, int n) {
    Eigen::MatrixXd u(1, n+1);
    for (int i = 0 ; i < n+1 ; i++) {
        u(0,i) = pow(-1.0, i)*pow(omega, 2*i);
    }
    return u;
}


bool does_ke_exist(ibex::IntervalVector& even_coeffs, Eigen::MatrixXd ve) {

    // std::cout << "prepare linprog 1" << std::endl;

    int n = ve.cols(); // no+1
    int m = ve.rows(); // ne+1
    Eigen::MatrixXd A(n+2*m, m);
    Eigen::VectorXd B(n+2*m   );
    Eigen::VectorXd c(1);
    Eigen::VectorXd f(m);

    //std::cout << "prepare linprog 2" << std::endl;
    for (int i = 0 ; i < m ; i ++) {
        B(i    ) =  0                  ;
        B(i+m  ) =  even_coeffs[i].ub();
        B(i+m*2) = -even_coeffs[i].lb();
        f(i    ) = 1;
    }
    
    //std::cout << "prepare linprog 3" << std::endl;

    A.block(0  , 0, n, m) =  ve;
    A.block(n  , 0, m, m) =  Eigen::MatrixXd::Identity(m,m);
    A.block(n+m, 0, m, m) = -Eigen::MatrixXd::Identity(m,m);
    int k = n+2*m;
    //std::cout << A << std::endl;
    //std::cout << B << std::endl;
    //std::cout << "start linprog " << std::endl;
    bool res = igl::linprog(c, A, B, k, f);
    //std::cout << "end linprog " << std::endl;
    //std::cout << "f: " << std::endl; std::cout << f << std::endl;
    return res;
}
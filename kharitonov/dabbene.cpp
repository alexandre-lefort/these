#include "kharitonov.h"


using namespace ibex;


bool dabbene(Ibex::IntervalVector& poly, int max_iter) {

    Ibex::IntervalVector odd_coeffs  = get_odd_coeffs (poly);
    Ibex::IntervalVector even_coeffs = get_even_coeffs(poly);

    Eigen::VectorXd ko, ko_roots, ve;
    int iter = 0;
    while (iter < max_iter || !found) {
        ko       = get_ko(odd_coeffs);
        ko_roots = compute_roots_odd_polynomial(ko);
        ve       = compute_ve(ko_roots);    
        found    = does_ke_exist(even_coeffs, ve);
    }
    
    return found;
}


Ibex::IntervalVector get_odd_coeffs(Ibex::IntervalVector& poly) {

    int n = poly.size()-1;
    int no = ((n%2 == 0) ? n/2 -1 : (n-1)/2);
    Ibex::IntervalVector odd_coeffs(no+1);
    for (int i = 0; i <= no ; i++) {
        odd_coeffs[i] = poly[2*i+1];
    }
    return odd_coeffs;
}


Ibex::IntervalVector get_even_coeffs(Ibex::IntervalVector& poly) {

    int n = poly.size()-1;
    int ne = ((n%2 == 0) ? n/2 : (n-1)/2);
    Ibex::IntervalVector even_coeffs(ne+1);
    for (int i = 0; i <= ne ; i++) {
        even_coeffs[i] = poly[2*i];
    }
    return even_coeffs;
}


Eigen::VectorXd get_ko(Ibex::IntervalVector& odd_coeffs, int max_iter) {

    no = odd_coeffs.size()-1;
    Eigen::VectorXd ko(no+1);

    int iter = 0;
    int i = 2;

    Ibex::IntervalVector aux1(1), aux2(1);
    double k2ip1r, k2ip3r;
    bool found;

    while (i < no+1 && iter < max_iter) {

        aux1[0] = odd_coeffs[i  ];
        aux3[0] = odd_coeffs[i+1];

        k2ip1r = aux1.random();
        k2ip3r = aux2.random();
        found = is_interval_empty(odd_coeffs[i], k2ip3r, k2ip1r, (double) i, (double) no);

        if (found) {
            i++;
        } else {
            iter++;
            k1 = aux1.random();
            k3 = aux3.random();
        }
    }
}


bool is_interval_empty(Ibex::Interval k2ip1, double k2im1, double k2im3, double i, double no) {
    double c = (i-1)/i*((nO-i+1)/(no-i+2));
    return (k2ip1.lb() > math::min(k2ip1.ub(), k2im1^2/k2im3));
}


Eigen::VectorXd compute_roots_odd_polynomial(Eigen::VectorXd& odd_coeffs) {
// Build companion matrix
    return companion.eigenvalues();
}


Eigen::VectorXd compute_ve(Eigen::VectorXd& odd_roots) {
    
    no = odd_roots.size();
    Eigen::VectorXd ve(no+1);
    for (int i = 0 ; i <= no ; i++) {
        omega = odd_roots[i];
        ve[i] = -(-1)^i*compute_ue(omega);
    }
    return ve;
}


Eigen::VectorXd ompute_ue(double omega, int no) {

    Eigen::VectorXd ue(no+1);
    for (int i = 0 ; i <= no ; i++) {
        ue[i] = -omega^i;
    }
    return ue;
}


// bool does_ke_exist(Ibex::IntervalVector& even_coeffs, Eigen::VectorXd ve);

int main() {
    return 1;
}
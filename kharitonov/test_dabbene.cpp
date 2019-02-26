#include "dabbene.h"


using namespace std;


void test_dabbene() {
    int max_iter = 1000;
    ibex::IntervalVector poly(5);
    poly[0] = ibex::Interval(1,2);
    poly[1] = ibex::Interval(3,4);
    poly[2] = ibex::Interval(5,6);
    poly[3] = ibex::Interval(7,8);
    poly[4] = ibex::Interval(9,10);
    bool res = dabbene(poly, max_iter);
    cout << res << endl;
}


void test_get_odd_coeffs() {
    ibex::IntervalVector poly(7);
    poly[0] = ibex::Interval(1,2);
    poly[1] = ibex::Interval(3,4);
    poly[2] = ibex::Interval(5,6);
    poly[3] = ibex::Interval(7,8);
    poly[4] = ibex::Interval(9,10);
    poly[5] = ibex::Interval(9,10);
    poly[6] = ibex::Interval(9,10);
    ibex::IntervalVector res = get_odd_coeffs(poly);
    cout << res << endl;
}


void test_get_even_coeffs() {
    ibex::IntervalVector poly(7);
    poly[0] = ibex::Interval(1,2);
    poly[1] = ibex::Interval(3,4);
    poly[2] = ibex::Interval(5,6);
    poly[3] = ibex::Interval(7,8);
    poly[4] = ibex::Interval(9,10);
    poly[5] = ibex::Interval(9,10);
    poly[6] = ibex::Interval(9,10);
    ibex::IntervalVector res = get_even_coeffs(poly);
    cout << res << endl;
}


void test_get_ko() {
    int max_iter = 1000;
    ibex::IntervalVector poly(7);
    poly[0] = ibex::Interval(1,2);
    poly[1] = ibex::Interval(3,4);
    poly[2] = ibex::Interval(5,6);
    poly[3] = ibex::Interval(7,8);
    poly[4] = ibex::Interval(-19,-10);
    poly[5] = ibex::Interval(-19,-10);
    poly[6] = ibex::Interval(9,10);
    ibex::IntervalVector aux = get_odd_coeffs(poly);
    Eigen::MatrixXd res = get_ko(aux, max_iter);
    cout << res << endl;
}


void test_compute_roots_odd_polynomial() {
    int max_iter = 1000;
    ibex::IntervalVector poly(7);
    poly[0] = ibex::Interval(1,2);
    poly[1] = ibex::Interval(3,4);
    poly[2] = ibex::Interval(5,6);
    poly[3] = ibex::Interval(7,8);
    poly[4] = ibex::Interval(-19,-10);
    poly[5] = ibex::Interval(-19,-10);
    poly[6] = ibex::Interval(9,10);
    ibex::IntervalVector aux = get_odd_coeffs(poly);
    Eigen::MatrixXd aux_2 = get_ko(aux, max_iter);
    Eigen::MatrixXd res = compute_roots_odd_polynomial(aux_2);
    cout << res << endl;
}


void test_compute_ve() {
    int max_iter = 1000;
    ibex::IntervalVector poly(7);
    poly[0] = ibex::Interval(1,2);
    poly[1] = ibex::Interval(3,4);
    poly[2] = ibex::Interval(5,6);
    poly[3] = ibex::Interval(7,8);
    poly[4] = ibex::Interval(-19,-10);
    poly[5] = ibex::Interval(-19,-10);
    poly[6] = ibex::Interval(9,10);
    ibex::IntervalVector aux = get_odd_coeffs(poly);
    ibex::IntervalVector auxe = get_even_coeffs(poly);
    Eigen::MatrixXd aux_2 = get_ko(aux, max_iter);
    Eigen::MatrixXd aux_3 = compute_roots_odd_polynomial(aux_2);
    Eigen::MatrixXd res = compute_ve(aux_3,auxe.size() -1);
    cout << res << endl;
}

void test_does_ke_exist() {
    int max_iter = 1000;
    ibex::IntervalVector poly(7);
    poly[0] = ibex::Interval(1,2);
    poly[1] = ibex::Interval(3,4);
    poly[2] = ibex::Interval(5,6);
    poly[3] = ibex::Interval(7,8);
    poly[4] = ibex::Interval(-19,-10);
    poly[5] = ibex::Interval(-19,-10);
    poly[6] = ibex::Interval(9,10);
    ibex::IntervalVector aux = get_odd_coeffs(poly);
    ibex::IntervalVector auxe = get_even_coeffs(poly);
    Eigen::MatrixXd aux_2 = get_ko(aux, max_iter);
    Eigen::MatrixXd aux_3 = compute_roots_odd_polynomial(aux_2);
    Eigen::MatrixXd aux_4 = compute_ve(aux_3,auxe.size()-1);
    bool res = does_ke_exist(auxe, aux_4);
    cout << res << endl;
}

int main() {
    test_dabbene();
    return 1;
}

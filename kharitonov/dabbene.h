#include "ibex.h"
#include "Eigen/Eigen"


bool                  dabbene(Ibex::IntervalVector& poly);
Ibex::IntervalVector  get_odd_coeffs (Ibex::IntervalVector& poly);
Ibex::IntervalVector  get_even_coeffs(Ibex::IntervalVector& poly);
Eigen::VectorXd       get_ko(Ibex::IntervalVector& odd_coeffs);
Ibex::Interval        is_interval_empty(double k2p1m, double k2p1p, double k2m1, double k2m3, int no);
Eigen::VectorXd       compute_roots_odd_polynomial(Eigen::VectorXd& odd_coeffs);
Eigen::VectorXd       compute_ve(Eigen::VectorXd& odd_roots);
bool                  does_ke_exist(Ibex::IntervalVector& even_coeffs, Eigen::VectorXd ve);
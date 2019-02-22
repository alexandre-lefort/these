#include "ibex.h"
#include "Eigen/Eigen"


Eigen::MatrixXd build_hurwitz_matrix_old(Eigen::VectorXd& poly);

int check_stability(Eigen::VectorXd& v);

int kharitonov(Ibex::IntervalVector& box);
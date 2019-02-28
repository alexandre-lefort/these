#ifndef KHARITONOV_H
#define KHARITONOV_H


#include "ibex.h"
#include <Eigen/Eigen>


Eigen::MatrixXd build_hurwitz_matrix_old(Eigen::VectorXd& poly);
bool check_stability(Eigen::VectorXd& v);
bool kharitonov(ibex::IntervalVector& box);


#endif
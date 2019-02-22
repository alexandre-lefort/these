#include "ibex.h"


using namespace ibex;

//std::pair<IntervalVector, IntervalVector> bissect_max(IntervalVector& inter) {
//     return (inter.bisect(std::max(inter.diam(),0.5));
//}

void test_kharitonov() {
    Variable x(5);
    IntervalVector x_ini(5);
    x_ini[0] = Interval(-10,10);
    x_ini[1] = Interval(-10,10);
    x_ini[2] = Interval(-10,10);
    x_ini[3] = Interval(-10000,10000);
    x_ini[4] = Interval(-10000,10000);

    Function coeffs1 = Function(x, x[3] - (x[0]*x[0] - 5*x[1] + 3)/(x[2]*0.25 + 25));
    Function coeffs2 = Function(x, 2*x[0]-5-x[1]);
    Function coeffs3 = Function(x, -x[3] + (x[0]*x[0] - 5*x[1] + 3)/(x[2]*0.25 + 25));
    Function coeffs4 = Function(x, -x[4] + (x[1]*x[0] - 5*x[2] + 3));

    SystemFactory fac;
    fac.add_var(x,x_ini);
    fac.add_ctr(coeffs2);
    //fac.add_ctr(coeffs2);
    //fac.add_ctr(coeffs3);
    //fac.add_ctr(coeffs4);

    NormalizedSystem* sys = new NormalizedSystem(fac);
   
    std::cout << "start" << std::endl;
    Set set(*sys,0.1);
    std::cout << "end" << std::endl;
    std::cout << set << std::endl;

}



int main() {
	test_kharitonov();
	return 1;
}

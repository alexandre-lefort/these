#include "kharitonov.h"


void test_kharitonov() {
    int max_iter = 1000;
    ibex::IntervalVector poly(5);
    poly[0] = ibex::Interval(1,2);
    poly[1] = ibex::Interval(3,4);
    poly[2] = ibex::Interval(5,6);
    poly[3] = ibex::Interval(7,8);
    poly[4] = ibex::Interval(9,10);
    bool res = kharitonov(poly);
    std::cout << res << std::endl;
}



int main() {
	test_kharitonov();
	return 1;
}

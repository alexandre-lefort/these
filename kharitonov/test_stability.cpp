#include "dabbene.h"
#include "kharitonov.h"

using namespace std;


void test_dabene() {
    int max_iter = 1000;
    ibex::IntervalVector poly(8);
    poly[0] = ibex::Interval(1,2);
    poly[1] = ibex::Interval(3,4);
    poly[2] = ibex::Interval(5,6);
    poly[3] = ibex::Interval(7,8);
    poly[4] = ibex::Interval(9,10);
    poly[5] = ibex::Interval(9,10);
    poly[6] = ibex::Interval(9,10);
    poly[7] = ibex::Interval(9,10);
    bool res = dabbene(poly, max_iter);
    cout << res << endl;
}

int main() {
    test_dabene();
    return 1;
}

#include "dabbene.h"
#include "kharitonov.h"


using namespace std;


void test_stability() {
    int iter_test = 1000;
    int max_iter = 1000;
    int degree   = 5;
    ibex::IntervalVector poly(degree+1);
    ibex::IntervalVector v(1);
    v[0] = ibex::Interval(0,10);
    
    for (int t = 0 ; t < iter_test ; t++)
    {
    for (int i = 0 ; i < degree+1 ; i++)
    {   
        double l1 = v.random()[0];
        double l2 = v.random()[0];
        if (l1 < l2) poly[i] = ibex::Interval(l1,l1);
        else         poly[i] = ibex::Interval(l1,l1);
    }

    bool res_dabbene = false;
    bool res_kharitonov = kharitonov(poly);
    if(res_kharitonov) {
    res_dabbene = dabbene(poly, max_iter);
    }
    if (res_kharitonov) cout << poly << " " << res_kharitonov << " " << res_dabbene << endl;
    }
}

int main() {
    test_stability();
    return 1;
}

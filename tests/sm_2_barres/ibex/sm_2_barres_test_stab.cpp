#include "ibex.h"
#include <string>
#include <ctime>
#include "ibex_Random.h"


using namespace std;
using namespace ibex;


void sm_2_barres() {

    Variable x(9);
    IntervalVector x_ini(9);
    
    x_ini[0] = Interval(-10,10);
    x_ini[1] = Interval(-10,10);
    x_ini[2] = Interval(-10,10);
    x_ini[3] = Interval(-10,10);
    x_ini[4] = Interval(-10,10);
    x_ini[5] = Interval(-10,10);
    x_ini[6] = Interval(-10,10);
    x_ini[7] = Interval(-10,10);
    x_ini[8] = Interval(-10,10);
    
    double CzW0 = -2.7  ;
    double CmQ0 = -0.38 ;
    double eps  =  0.0000000001  ;
    
    Variable y(3);
    IntervalVector y_ini(3);
    
    y_ini[0] = Interval(CzW0 - eps*fabs(CzW0), CzW0 + eps*fabs(CzW0));
    y_ini[1] = Interval(CmQ0 - eps*fabs(CmQ0), CmQ0 + eps*fabs(CmQ0));
    y_ini[2] = Interval(-3,1);    

    Variable p(2);
    IntervalVector p_ini(2);
    
    p_ini[0] = Interval(CzW0 - eps*fabs(CzW0), CzW0 + eps*fabs(CzW0));
    p_ini[1] = Interval(CmQ0 - eps*fabs(CmQ0), CmQ0 + eps*fabs(CmQ0));

    Function crit_hinf_zz1("../functions/Tzz1.txt" )       ; cout << "hinf criteria loaded" << endl;
    Function crit_hinf_zz2("../functions/Tzz2.txt" )       ; cout << "hinf criteria loaded" << endl;
    Function crit_hinf_zb1("../functions/Tzb1.txt" )       ; cout << "hinf criteria loaded" << endl;
    Function crit_hinf_zb2("../functions/Tzb2.txt" )       ; cout << "hinf criteria loaded" << endl;
    
    Function coeffs("../functions/Tstab_coefs.txt");

    const ExprNode& stab1 = coeffs(x,p)[0];
    const ExprNode& stab2 = coeffs(x,p)[1];
    const ExprNode& stab3 = coeffs(x,p)[2];
    const ExprNode& stab4 = coeffs(x,p)[3];
    const ExprNode& stab5 = coeffs(x,p)[4];
    const ExprNode& stab6 = coeffs(x,p)[5];
    
    Function stab(x,p,ibex::max(ibex::max(ibex::max(ibex::max(-stab6,-stab4),-stab2),ibex::pow(stab1,2)*ibex::pow(stab6,2) + ibex::pow(stab2,2)*ibex::pow(stab5,2) + stab1*ibex::pow(stab4,2)*stab5 + stab2*ibex::pow(stab3,2)*stab6 - 2*stab1*stab2*stab5*stab6 - stab1*stab3*stab4*stab6 - stab2*stab3*stab4*stab5),stab1*stab4 - stab2*stab3));

    Function f_max1(x,y,crit_hinf_zz1(x,y));
    Function f_max2(x,y,crit_hinf_zz2(x,y));
    Function f_max3(x,y,crit_hinf_zb1(x,y));
    Function f_max4(x,y,crit_hinf_zb2(x,y));
        
    Function goals(x,y,ibex::max(f_max1(x,y),(ibex::max(f_max2(x,y),ibex::max(f_max3(x,y),f_max4(x,y))))));
    
    IntervalVector inter(1);
    inter[0] = Interval(-10,9.9);

    for (int i = 0 ; i < 9 ; i++) {

        IntervalVector v_test(12);

        for (int j = 0 ; j < 9 ; j++) {
            double v1 = (double) i;
            double v2 = v1 + 1.0;
            v_test[j] = Interval(v1,v2);
        }

        IntervalVector inter2(1);
        inter2[0] = Interval(-3,0.999);
        double v3 = i;
        v_test[ 9] = y_ini[0];
        v_test[10] = y_ini[1];
        v_test[11] = Interval(v3 + 0.1);

        cout << stab.eval(v_test)  << endl;
        //cout << goals.eval(v_test) << endl;
    }
    

}


int main(int argc, char* argv[]) {

RNG::srand(1234);
sm_2_barres();

}

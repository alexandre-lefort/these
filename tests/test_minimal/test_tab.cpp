#include "ibex.h"
#include <string>
#include <ctime>


using namespace std;
using namespace ibex;


void test_tab() {

    int num_thread = 16;

    // Problem formulation
    Variable x(2),y(2);

    IntervalVector x_ini(2);
    x_ini[0] = Interval(-100,100);
    x_ini[1] = Interval(-100,100);

    IntervalVector y_ini(2);
    y_ini[0] = Interval(-5,5);
    y_ini[1] = Interval(-5,5);

    double x_prec(1e-12),y_prec(1e-10),stop_prec(0.1);

    Function goal = Function(x,y,4*ibex::pow(x[0] - 2,2) - 2*ibex::pow(y[0],2) + ibex::pow(x[0],2)*y[0] - ibex::pow(y[1],2) + 2*ibex::pow(x[1],2)*y[1]);

    SystemFactory fac_x;
    fac_x.add_var(x,x_ini);

    Function stab = Function(x,y, x[0] + y[0] - 1000);

    cout<<"constraints loaded"<<endl;

    NormalizedSystem sys_x_1(fac_x);

    SystemFactory fac_xy;
    fac_xy.add_var(x,x_ini);
    fac_xy.add_var(y,y_ini);
    fac_xy.add_goal(goal);
    cout<<"systems ok"<<endl;
    NormalizedSystem sys_xy_1(fac_xy);

    SystemFactory fac_fa_y_sys;
    fac_fa_y_sys.add_var(x,x_ini);
    fac_fa_y_sys.add_var(y,y_ini);
    fac_fa_y_sys.add_goal(stab);

    double prec_fa_y = 1e-3;

    NormalizedSystem fa_y_sys_1(fac_fa_y_sys);

    vector<System>      sys_x    = vector<System>(); 
    vector<System>      sys_xy   = vector<System>(); 
    vector<CtcIdentity> xy_ctc   = vector<CtcIdentity>();
    vector<CtcIdentity> x_ctc_id = vector<CtcIdentity>(); 
    vector<CtcIdentity> fa_y_ctc = vector<CtcIdentity>(); 

    for (int i = 0; i<num_thread ; i++) {

        sys_x.push_back(System(sys_x_1));
        sys_xy.push_back(System(sys_xy_1));
        xy_ctc.push_back(CtcIdentity(x_ini.size() + y_ini.size()));
        x_ctc_id.push_back(CtcIdentity(x_ini.size()));
        fa_y_ctc.push_back(CtcIdentity(x_ini.size()+y_ini.size()));
    }

    for (int i = 0 ; i < num_thread ; i++) {
        //std::cout << sys_x[i].goal   << std::endl;
        //std::cout << sys_xy[i].goal   << std::endl;
        std::cout << xy_ctc[i].nb_var   << std::endl;
        std::cout << x_ctc_id[i].nb_var << std::endl;
        std::cout << fa_y_ctc[i].nb_var << std::endl;
    }
}



int main(int argc, char* argv[]) {
    test_tab();
}

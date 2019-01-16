#include "ibex.h"
#include <string>

using namespace std;
using namespace ibex;

void ex_robsyn_phd(){
    Variable x(2),y(3);
    IntervalVector x_ini(2,Interval(0.1,5));
    x_ini[0] = Interval(0.1,5);
    x_ini[1] = Interval(0.1,5);
    IntervalVector y_ini(3,Interval(-3,3));
//    x_ini[3] = Interval(0,5);
    y_ini[0] = Interval(0.5,1);
    y_ini[1] = Interval(0.5,1);
    y_ini[2] = Interval(-2,2);
    double x_prec(1e-12),y_prec(1e-10),stop_prec(0.1);

    Function goal("objective.txt");

    cout<<"functions loaded"<<endl;

    SystemFactory fac_x;
    fac_x.add_var(x,x_ini);
    CtcIdentity x_ctc_id(x_ini.size());

    Variable p(2);
    IntervalVector p_ini(2);
    p_ini[0] = y_ini[0];
    p_ini[1] = y_ini[1];
    //fac_x.add_var(p,p_ini);

    Function stab("csts.txt");


    Function c1f(x,p,-stab(x,p)[0]);
//    Function c2f(x,p,-stab(x,p)[1]);
//    Function c3f(x,p,-stab(x,p)[2]);

    cout<<"constraints loaded"<<endl;

    NormalizedSystem sys_x(fac_x);

    SystemFactory fac_xy;
    fac_xy.add_var(x,x_ini);
    fac_xy.add_var(y,y_ini);
    fac_xy.add_goal(goal);
    CtcIdentity xy_ctc(x_ini.size() + y_ini.size());
    cout<<"systems ok"<<endl;
    NormalizedSystem sys_xy(fac_xy);

    SystemFactory fac_fa_y_sys;
    fac_fa_y_sys.add_var(x,x_ini);
    fac_fa_y_sys.add_var(p,p_ini);
//    Function max12(x,p,ibex::max(c1f(x,p),c2f(x,p)));
//    Function max3(x,p,ibex::max(max12(x,p),c3f(x,p)));
    fac_fa_y_sys.add_goal(c1f);

    double prec_fa_y = 1e-3;

    NormalizedSystem fa_y_sys(fac_fa_y_sys);
    CtcIdentity fa_y_ctc(x_ini.size()+p_ini.size());


    OptimMinMax oo(sys_x, sys_xy, fa_y_sys, x_ctc_id, xy_ctc, fa_y_ctc, x_prec, y_prec,stop_prec,prec_fa_y);
//    OptimMinMax oo(sys_x, sys_xy,x_ctc_id,xy_ctc,x_prec, y_prec,stop_prec);
    oo.list_elem_absolute_max = 0;
    oo.list_rate = 0;
    oo.critpr = 0.4;
    oo.heap_prob = 0.1;
    oo.min_prec_coef = 0;
    oo.iter = 15;
    oo.local_iter = 0;
    oo.visit_all = true;
    oo.nb_point = 1;

    oo.list_elem_absolute_max_csp = 0;
    oo.iter_csp = 10;
    oo.critpr_csp = 0.3;
    oo.list_rate_csp = 0;
    oo.min_prec_coef_csp = 0;
    oo.local_iter_csp = 0;
    oo.visit_all_csp = true;

    oo.trace=1;
    oo.trace_freq = 1;
    oo.timeout=3600;
    Optim::Status res = oo.optimize(x_ini);
    oo.report();

}





int main(int argc, char* argv[]) {

ex_robsyn_phd();

}

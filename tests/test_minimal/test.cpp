#include "ibex.h"
#include <string>
#include <ctime>
#include "ibex_Random.h"


using namespace std;
using namespace ibex;


void problem_2() {

    Variable x(2),y(2);

    IntervalVector x_ini(2);
    x_ini[0] = Interval(-100,100);
    x_ini[1] = Interval(-100,100);

    IntervalVector y_ini(2);
    y_ini[0] = Interval(-5,5);
    y_ini[1] = Interval(-5,5);

    int num_thread = 8;

    double x_prec(1e-6),y_prec(1e-6),stop_prec(0.001);

    std::vector<Function*>         goals    = std::vector<Function*>();   
    std::vector<Function*>         stabs    = std::vector<Function*>();   
    std::vector<NormalizedSystem*> sys_x    = std::vector<NormalizedSystem*>();
    std::vector<NormalizedSystem*> sys_xy   = std::vector<NormalizedSystem*>();
    std::vector<NormalizedSystem*> fa_y_sys = std::vector<NormalizedSystem*>();
    std::vector<Ctc*>              x_ctc_id = std::vector<Ctc*>();          
    std::vector<Ctc*>              xy_ctc   = std::vector<Ctc*>();          
    std::vector<Ctc*>              fa_y_ctc = std::vector<Ctc*>();          
    std::cout << "fill structures" << std::endl;

    for (int i = 0 ; i <num_thread ; i++) {

        goals.push_back(new Function(x,y,4*ibex::pow(x[0] - 2,2) - 2*ibex::pow(y[0],2) + ibex::pow(x[0],2)*y[0] - ibex::pow(y[1],2) + 2*ibex::pow(x[1],2)*y[1]));

        SystemFactory fac_x;
        fac_x.add_var(x,x_ini);

        stabs.push_back(new Function(x,y, x[0] + y[0] - 1000));
        
        SystemFactory fac_xy;
        fac_xy.add_var(x,x_ini);
        fac_xy.add_var(y,y_ini);
        fac_xy.add_goal(*goals[i]);
       
        SystemFactory fac_fa_y_sys;
        fac_fa_y_sys.add_var(x,x_ini);
        fac_fa_y_sys.add_var(y,y_ini);
        fac_fa_y_sys.add_goal(*stabs[i]);


        std::cout << i << std::endl;
        sys_x.push_back(     new NormalizedSystem(fac_x));
        sys_xy.push_back(    new NormalizedSystem(fac_xy));
        fa_y_sys.push_back(  new NormalizedSystem(fac_fa_y_sys));
        std::cout << "sys ok" << sys_x.size() << std::endl;
        xy_ctc.push_back(    new CtcHC4(*sys_xy[i]));
        x_ctc_id.push_back(  new CtcIdentity(x_ini.size()));
        fa_y_ctc.push_back(  new CtcIdentity(x_ini.size()+y_ini.size()));
        std::cout << "ctc ok" << xy_ctc.size() << std::endl;
    }

    double prec_fa_y = 1e-3;

    OptimMinMax oo(sys_x, sys_xy, fa_y_sys, x_ctc_id, xy_ctc, fa_y_ctc, x_prec, y_prec,stop_prec,prec_fa_y, num_thread);
    

    oo.list_elem_absolute_max = 0;
    oo.list_rate = 0;
    oo.critpr = 0.4;
    oo.heap_prob = 0.1;
    oo.min_prec_coef = 0;
    oo.iter = 100;
    oo.local_iter = 0;
    oo.visit_all = false;
    oo.nb_point = 1;

    oo.list_elem_absolute_max_csp = 0;
    oo.iter_csp = 100;
    oo.critpr_csp = 0.3;
    oo.list_rate_csp = 0;
    oo.min_prec_coef_csp = 0;
    oo.local_iter_csp = 0;
    oo.visit_all_csp = false;

    oo.trace=1;
    oo.trace_freq = 1;
    oo.timeout=3600;


    Optim::Status res = oo.optimize(x_ini);

    oo.report();

}

void problem_1() {

    Variable x(1),y(1);

    IntervalVector x_ini(1);
    x_ini[0] = Interval(-100,100);

    IntervalVector y_ini(1);
    y_ini[0] = Interval(-5,5);

    double x_prec(1e-12),y_prec(1e-10),stop_prec(0.001);

    int num_thread = 2;

    std::vector<Function*>         goals    = std::vector<Function*>();   
    std::vector<Function*>         stabs    = std::vector<Function*>();   
    std::vector<NormalizedSystem*> sys_x    = std::vector<NormalizedSystem*>();
    std::vector<NormalizedSystem*> sys_xy   = std::vector<NormalizedSystem*>();
    std::vector<NormalizedSystem*> fa_y_sys = std::vector<NormalizedSystem*>();
    std::vector<Ctc*>              x_ctc_id = std::vector<Ctc*>();          
    std::vector<Ctc*>              xy_ctc   = std::vector<Ctc*>();          
    std::vector<Ctc*>              fa_y_ctc = std::vector<Ctc*>();          

    std::cout << "fill structures" << std::endl;

    for (int i = 0 ; i< num_thread ; i++) {

        goals.push_back(new Function(x,y,ibex::pow(x[0] - 2,2) + ibex::pow(x[0],2)*y[0]));

        SystemFactory fac_x;
        fac_x.add_var(x,x_ini);

        stabs.push_back(new Function(x,y, x[0] + y[0] - 1000));

        SystemFactory fac_xy;
        fac_xy.add_var(x,x_ini);
        fac_xy.add_var(y,y_ini);
        fac_xy.add_goal(*goals[i]);

        SystemFactory fac_fa_y_sys;
        fac_fa_y_sys.add_var(x,x_ini);
        fac_fa_y_sys.add_var(y,y_ini);
        fac_fa_y_sys.add_goal(*stabs[i]);

        std::cout << i << std::endl;
        sys_x.push_back(     new NormalizedSystem(fac_x));
        sys_xy.push_back(    new NormalizedSystem(fac_xy));
        fa_y_sys.push_back(  new NormalizedSystem(fac_fa_y_sys));
        std::cout << "sys ok1" << sys_x.size() << std::endl;
        xy_ctc.push_back(    new CtcHC4(*sys_xy[i]));

        cout << *fa_y_sys[i] << endl;

        std::cout << "sys ok2" << sys_x.size() << std::endl;
        x_ctc_id.push_back(  new CtcIdentity(x_ini.size()));
        std::cout << "sys ok3" << sys_x.size() << std::endl;
        fa_y_ctc.push_back(  new CtcIdentity(x_ini.size() + y_ini.size()));
        std::cout << "ctc ok" << xy_ctc.size() << std::endl;
    }

    
    double prec_fa_y = 1e-3;
    OptimMinMax oo(sys_x, sys_xy, fa_y_sys, x_ctc_id, xy_ctc, fa_y_ctc, x_prec, y_prec, stop_prec, prec_fa_y, num_thread);
    std::cout << "optim object created" << std::endl;
    //oo.num_thread = 4;
    oo.list_elem_absolute_max = 0;
    oo.list_rate = 0;
    oo.critpr = 0.4;
    oo.heap_prob = 0.1;
    oo.min_prec_coef = 1e-3;
    oo.iter = 100;
    oo.local_iter = 0;
    oo.visit_all = false;
    oo.nb_point = 1;
    
    oo.list_elem_absolute_max_csp = 0;
    oo.iter_csp = 100;
    oo.critpr_csp = 0.3;
    oo.list_rate_csp = 0;
    oo.min_prec_coef_csp = 1e-3;
    oo.local_iter_csp = 0;
    oo.visit_all_csp = false;
    
    oo.trace=1;
    oo.trace_freq = 1;
    oo.timeout=3600;
    
    std::clock_t start;
    double duration;
    
    start = std::clock();
    
    std::cout << "optim start" << std::endl;
    Optim::Status res = oo.optimize(x_ini);
    
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    
    std::cout<<"printf: "<< duration <<'\n';
    
    oo.report();
   

}



int main(int argc, char* argv[]) {

RNG::srand(1234);
problem_2();

}

//============================================================================
//                                  I B E X
// File        : ibex_OptimMinMax.cpp
// Author      : Dominique Monnet, Jordan Ninin
// License     : See the LICENSE file
// Created     : Oct 1, 2016
//============================================================================

#include "ibex_OptimMinMax.h"
#include "ibex_LargestFirst.h"
#include "ibex_Timer.h"
#include <stdio.h>
#include "ibex_DataMinMax.h"
#include "ibex_NoBisectableVariableException.h"
#include "ibex_SystemFactory.h"

#define NOT_SATISFIED  0
#define NOT_DETERMINED 1
#define SATISFIED      2


using namespace std;

namespace ibex {

//********* Default parameters for light optim minmax solver ************
const int OptimMinMax::default_iter                   = 10    ;
const int OptimMinMax::default_list_rate              = 0     ;
const double OptimMinMax::default_min_prec_coef       = 10    ;
const int OptimMinMax::default_list_elem_absolute_max = 500   ;
const int OptimMinMax::default_prob_heap              = 10    ; //10% to pop second heap in light_solver
const bool OptimMinMax::default_visit_all             = false ;

const int OptimMinMax::default_nb_point       = 1  ;
const double OptimMinMax::default_perf_thresh = 0.3;

// Csp default parameters for light solver
const int OptimMinMax::default_iter_csp                   = 10;
const int OptimMinMax::default_list_rate_csp              = 0;
const double OptimMinMax::default_min_prec_coef_csp       = 100;
const int OptimMinMax::default_list_elem_absolute_max_csp = 100;
const int OptimMinMax::default_prob_heap_csp              = 0; //0% to pop second heap in light_solver
const bool OptimMinMax::default_visit_all_csp             = false;



OptimMinMax::OptimMinMax(std::vector<NormalizedSystem*>  x_sys_t        , 
                         std::vector<NormalizedSystem*>  xy_sys_t       , 
                         std::vector<Ctc*>               x_ctc_t        , 
                         std::vector<Ctc*>               xy_ctc_t       ,
                         double                          prec_x         , 
                         double                          prec_y         , 
                         double                          goal_rel_prec  ,
                         int                             num_thread_t   ):
    Optim(x_sys_t[0]->nb_var, new CellDoubleHeap(*new CellCostFmaxlb_opt(), *new CellCostFmaxub_opt()),
          prec_x, goal_rel_prec, 0, 1), // attention meme precision en relatif et en absolue
    x_box_init(x_sys_t[0]->box),
    y_box_init(xy_sys_t[0]->box.subvector(x_sys_t[0]->nb_var, xy_sys_t[0]->nb_var-1)),
    y_box_init_fa(IntervalVector(1)),  
    propag(true),
    trace_freq(10000),
    num_thread(num_thread_t),
    list_rate(default_list_rate),
    list_elem_absolute_max(default_list_elem_absolute_max),
    iter(default_iter),
    min_prec_coef(default_min_prec_coef),
    critpr(default_prob_heap),
    list_rate_csp(0),                       // unused                   
    list_elem_absolute_max_csp(0),          // unused                  
    iter_csp(0),                            // unused
    min_prec_coef_csp(0.0),                 // unused           
    critpr_csp(0),                          // unused       
    only_csp(false),
    heap_prob(0),                           // unused
    visit_all(default_visit_all),
    visit_all_csp(default_visit_all_csp),   // unused
    nb_point(default_nb_point),
    perf_thresh(default_perf_thresh),
    prec_y(prec_y),
    fa_y_cst(false),
    min_goal(x_sys_t[0]->goal != NULL),
    prec_fa_y(0),
    x_ctc(x_ctc_t),                   
    x_sys(x_sys_t),      
    lsolve(std::vector<LightOptimMinMax*>()),            
    bsc(std::vector<Bsc*>()),                     
    minus_goal_y_at_x(std::vector<Function*>()),                    
    fa_lsolve(std::vector<LightOptimMinMax*>()),                
    minus_goal_csp_y_at_x(std::vector<Function*>())                  
{
    // TODO : test structures ExprNode in a unit test
    for (int i = 0 ; i < num_thread ; i ++) {
        lsolve.push_back(new LightOptimMinMax(*xy_sys_t[i], *xy_ctc_t[i]));
        bsc.push_back(new LargestFirst());  
        fa_lsolve.push_back(new LightOptimMinMax(*xy_sys_t[i], *xy_ctc_t[i], true)); // useless if no fa cst but need to construct it...
        
       if(!min_goal && xy_sys_t[0]->goal != NULL) {
   //        // goal function reformulation as min instead of max for local solver
           Array<const ExprNode> args((xy_sys_t[0]->goal)->nb_arg());
           Array<const ExprSymbol> var;
           for(int j = 0 ; j < (xy_sys_t[0]->goal)->nb_arg() ; j++) {
               const ExprSymbol& a = ExprSymbol::new_((xy_sys_t[i]->goal)->arg(j).dim);
               var.add(a);
               args.set_ref(j,a);
           }
           minus_goal_y_at_x.push_back(new Function(var,-(*(xy_sys_t[i]->goal))(args)));    
           Affine2Eval* aff_eval = new Affine2Eval(*(xy_sys_t[i]->goal));
           lsolve[i]->affine_goal = aff_eval;
           lsolve[i]->goal_abs_prec = goal_rel_prec/100; // set goal prec of maximization problem lower than minimization
           fa_lsolve[i]->goal_abs_prec = 1e-2;
       }
   }

   omp_init_lock(&bufferlock); // omp
};

 // Todo :
OptimMinMax::OptimMinMax(std::vector<NormalizedSystem*> x_sys_t            , 
                         std::vector<NormalizedSystem*> xy_sys_t           , 
                         std::vector<NormalizedSystem*> max_fa_y_cst_sys_t , 
                         std::vector<Ctc*>              x_ctc_t            , 
                         std::vector<Ctc*>              xy_ctc_t           , 
                         std::vector<Ctc*>              y_fa_ctc_t         ,
                         double                         prec_x             , 
                         double                         prec_y             , 
                         double                         goal_rel_prec      , 
                         double                         fa_cst_prec        ,
                         int                            num_thread_t       ):
    Optim(x_sys_t[0]->nb_var, new CellDoubleHeap(*new CellCostFmaxlb_opt(), *new CellCostFmaxub_opt()),
          prec_x, goal_rel_prec, goal_rel_prec, 1), // attention meme precision en relatif et en absolue
    x_box_init(x_sys_t[0]->box),
    y_box_init(xy_sys_t[0]->box.subvector(x_sys_t[0]->nb_var, xy_sys_t[0]->nb_var-1)),
    y_box_init_fa(max_fa_y_cst_sys_t[0]->box.subvector(x_sys_t[0]->nb_var, max_fa_y_cst_sys_t[0]->nb_var-1)),
    propag(true),
    trace_freq(10000),
    num_thread(num_thread_t),
    list_rate(default_list_rate),
    list_elem_absolute_max(default_list_elem_absolute_max),
    iter(default_iter),
    min_prec_coef(default_min_prec_coef),
    critpr(default_prob_heap), 
    list_rate_csp(default_list_rate_csp),
    list_elem_absolute_max_csp(default_list_elem_absolute_max_csp),
    iter_csp(default_iter_csp),
    min_prec_coef_csp(default_min_prec_coef_csp),
    critpr_csp(default_prob_heap_csp),
    only_csp(false),
    heap_prob(0),
    visit_all(default_visit_all),
    visit_all_csp(default_visit_all_csp),
    nb_point(default_nb_point),
    perf_thresh(default_perf_thresh),
    prec_y(prec_y),
    fa_y_cst(true),
    min_goal(x_sys_t[0]->goal != NULL),
    prec_fa_y(fa_cst_prec),
    x_ctc(x_ctc_t),                   
    x_sys(x_sys_t),      
    lsolve(std::vector<LightOptimMinMax*>()),           
    bsc(std::vector<Bsc*>()),                     
    minus_goal_y_at_x(std::vector<Function*>()),                  
    fa_lsolve(std::vector<LightOptimMinMax*>()),                
    minus_goal_csp_y_at_x(std::vector<Function*>())            
{      
    for (int i = 0 ; i < num_thread ; i ++) {
        lsolve.push_back(new LightOptimMinMax(*xy_sys_t[i], *xy_ctc_t[i]));
        bsc.push_back(new LargestFirst()); 

        fa_lsolve.push_back(new LightOptimMinMax(*max_fa_y_cst_sys_t[i], *y_fa_ctc_t[i], true));
        
        if(!min_goal && xy_sys_t[i]->goal !=NULL) {
            Array<const ExprNode> args((xy_sys_t[i]->goal)->nb_arg());
            Array<const ExprSymbol> var;
            for(int j = 0 ; j < (xy_sys_t[i]->goal)->nb_arg() ; j++) {
                const ExprSymbol& a = ExprSymbol::new_((xy_sys_t[i]->goal)->arg(j).dim);
                var.add(a);
                args.set_ref(j,a);
            }
            minus_goal_y_at_x.push_back(new Function(var,-(*(xy_sys_t[i]->goal))(args)));    
            Affine2Eval* aff_eval = new Affine2Eval(*(xy_sys_t[i]->goal));
            lsolve[i]->affine_goal = aff_eval;
            lsolve[i]->goal_abs_prec = goal_rel_prec/100; // set goal prec of maximization problem lower than minimization
            fa_lsolve[i]->goal_abs_prec = 1e-2;
        }
     
        Array<const ExprNode> args_csp((xy_sys_t[i]->goal)->nb_arg());
        Array<const ExprSymbol> var_csp;
        for(int j = 0; j < (xy_sys_t[0]->goal)->nb_arg();j++) {
            const ExprSymbol& a = ExprSymbol::new_((xy_sys_t[i]->goal)->arg(j).dim);
            var_csp.add(a);
            args_csp.set_ref(j,a);
        }
     
        minus_goal_csp_y_at_x.push_back(new Function(var_csp,-(*(max_fa_y_cst_sys_t[i]->goal))(args_csp)));         
        Affine2Eval* aff_eval_csp = new Affine2Eval(*(max_fa_y_cst_sys_t[i]->goal));
        fa_lsolve[i]->affine_goal = aff_eval_csp;
        fa_lsolve[i]->goal_abs_prec = 1e-2;
    }
    omp_init_lock(&bufferlock); // omp
};

//Todo :
OptimMinMax::~OptimMinMax() {
    cout<<"Call destructor for optiminmax object: "<<this<<endl;
    cout<<"try to flush buffer"<<endl;
    if(buffer!=NULL) {
        buffer->flush();
        cout<<"buffer flushed"<<endl;
        delete &(buffer->cost1());
        delete &(buffer->cost2());
        cout<<"buffer cost func deleted"<<endl;
        delete buffer;
    }
    cout<<"buffer deleted"<<endl;

    omp_destroy_lock(&bufferlock); // omp

    for (int i = 0 ; i < num_thread ; i++) {
        delete(lsolve[i]);
        delete(bsc[i]);
        delete(minus_goal_y_at_x[i]);
        delete(fa_lsolve[i]); 
        delete(minus_goal_csp_y_at_x[i]); 
    }  
}


bool OptimMinMax::check_optimizer() {
    if(x_sys[0]->goal != NULL && lsolve[0]->xy_sys.goal != NULL) { // two objective functions-> error
        cout<<" Error: Two ojective functions found, choose either x depending function to minimize or x and y depending function for min max optim."<<endl;
        return false;
    }
    return true;
}


Optim::Status OptimMinMax::optimize() {
    bool problem_ok = check_optimizer();
    if(!problem_ok) { return INFEASIBLE; }
    return optimize(x_sys[0]->box,POS_INFINITY);
}


void OptimMinMax::init_lsolve() {
    for (int i = 0 ; i < num_thread ; i ++) {
        lsolve[i]->visit_all = visit_all;
    }
}


void OptimMinMax::init_fa_lsolve() {
    for (int i = 0 ; i < num_thread ; i ++) {
        fa_lsolve[i]->visit_all = visit_all_csp;
    }
}


Optim::Status OptimMinMax::optimize(const IntervalVector& x_box_ini1, double obj_init_bound) {
    
    cout<<"start optimization"<<endl;

    loup = obj_init_bound;
    uplo = NEG_INFINITY;
    uplo_of_epsboxes = POS_INFINITY;

    x_box_init = x_box_ini1;

    //****** initialization of the first Cell********
    Cell * root = new Cell(x_box_init);
    buffer->cost1().add_backtrackable(*root);
    buffer->cost2().add_backtrackable(*root);

    bsc[0]->add_backtrackable(*root);
    lsolve[0]->add_backtrackable(*root, y_box_init, critpr);
    if(fa_y_cst) {
        fa_lsolve[0]->add_backtrackable(*root, y_box_init_fa, critpr_csp);
    }

    buffer->critpr = heap_prob;

    init_lsolve()   ; //************* set light optim minmax solver param     *********   
    init_fa_lsolve(); //************* set csp light optim minmax solver param *********

    //****** x_heap initialization ********

    nb_cells = 0;
    buffer->flush();

    // *** initialisation Algo ***
    loup_changed = false;
    double ymax;
    initial_loup=obj_init_bound;
    loup_point = x_box_init.mid();
    time = 0;

    //monitoring variables, used to track upper bound, lower bound, number of elem in y_heap and heap_save at each iteration
    std::vector<double> ub, lb, nbel, nbyel;
    long long int nbel_count(0);
    vector<bool> handle_res = vector<bool>(num_thread);

    if(!handle_cell(root,0)) { return INFEASIBLE; }
    if (loup_changed) {
        ymax = compute_ymax();
        buffer->contract(ymax);

        if (ymax <= NEG_INFINITY) {
            if (trace) cout << " infinite value for the minimum " << endl;
        }
        if (trace) cout <<  "iter="<< nb_cells <<",  size_heap="<< buffer->size()<< ",  ymax=" << ymax << ",  uplo= " <<  uplo<< endl;
    }
    if (trace) cout <<  "iter="<< nb_cells <<",  size_heap="<< buffer->size()<< ",  ymax=" << ymax << ",  uplo= " <<  uplo<< endl;
    update_uplo();

    Timer::reset_time();
    Timer::start();

    //cout<<"start optimization loop"<<endl;

    try {
        while (!buffer->empty() && (loup-uplo)>goal_rel_prec) {

            if (trace >= 2) buffer->print(cout);
            loup_changed = false;
            Cell *c = buffer->pop();

            //====== only for purtpose of non inheritance comparison
            if(!propag){
                DataMinMax * data_x = &(c->get<DataMinMaxOpti>());

                data_x->y_heap->flush();
                Cell * y_cell = new Cell(y_box_init);
                y_cell->add<OptimData>();
                data_x->y_heap->push(y_cell);
                if(fa_y_cst) {
                    DataMinMax * data_x_csp = &(c->get<DataMinMaxCsp>());
                    data_x_csp->y_heap->flush();
                    Cell * y_cell_csp = new Cell(y_box_init_fa);
                    y_cell_csp->add<OptimData>();
                    data_x_csp->y_heap->push(y_cell_csp);
                }
            }

            try {
                // omp
                std::vector<Cell*> new_cells = nsect_cell(num_thread, c);
                delete c;
                
                #pragma omp parallel num_threads(num_thread)
                {
                    int ID = omp_get_thread_num();
                    handle_res[ID] = handle_cell(new_cells[ID], ID);
                }
                // end omp

                if (uplo_of_epsboxes == NEG_INFINITY) {
                    cout << " possible infinite minimum " << endl;
                    break;
                }
                
                if (loup_changed) {
                    ymax = compute_ymax();
                    buffer->contract(ymax);

                    if (ymax <= NEG_INFINITY) {
                        if (trace) cout << " infinite value for the minimum " << endl;
                        break;
                    }
                    if (trace) cout <<  "iter="<< nb_cells <<",  size_heap="<< buffer->size()<< ",  ymax=" << ymax << ",  uplo= " <<  uplo<< endl;
                }
                if (trace) cout <<  "iter="<< nb_cells <<",  size_heap="<< buffer->size()<< ",  ymax=" << ymax << ",  uplo= " <<  uplo<< endl;
                update_uplo();

                Timer::check(timeout);

            }
            catch (NoBisectableVariableException& ) {
                bool res = handle_cell(c,0);
                if (res) update_uplo_of_epsboxes(c->get<DataMinMaxOpti>().fmax.lb());
                //if (trace>=1) cout << "epsilon-box found: uplo cannot exceed " << uplo_of_epsboxes << endl;
                update_uplo(); // the heap has changed -> recalculate the uplo
            }
        }
    }
    catch (TimeOutException& ) {
        Timer::stop();
        time = Timer::get_time();
        return TIME_OUT;
    }

    Timer::stop();
    time = Timer::get_time();

    if (uplo_of_epsboxes == POS_INFINITY && (loup==POS_INFINITY || (loup==initial_loup && goal_abs_prec==0 && goal_rel_prec==0)))
        return INFEASIBLE;
    else if (loup==initial_loup)
        return NO_FEASIBLE_FOUND;
    else if (uplo_of_epsboxes == NEG_INFINITY)
        return UNBOUNDED_OBJ;
    else
        return SUCCESS;
}


std::vector<Cell*> OptimMinMax::nsect_cell(int n, Cell *c) {

    int nstep = std::floor(log2(n));

    std::vector<Cell*> cells_1 = std::vector<Cell*>();
    std::vector<Cell*> cells_2 = std::vector<Cell*>();

    cells_1.push_back(c);

    for (int i = 0 ; i < nstep ; i++) {
        for (int j = 0 ; j < std::pow(2,i) ; j ++) {
            pair<Cell*,Cell*> new_cells = bsc[0]->bisect_cell(*cells_1[j]);
            cells_2.push_back(new_cells.first);
            cells_2.push_back(new_cells.second);
        }
        cells_1 = cells_2;
        cells_2.clear();
    }

    return cells_1;
}


bool  OptimMinMax::handle_cell(Cell * x_cell, int ith) {

    DataMinMaxOpti * data_x = &(x_cell->get<DataMinMaxOpti>());

    //***************** contraction w.r.t constraint on x ****************
    IntervalVector tmpbox(x_cell->box);
    int res_cst = check_constraints(x_cell, false, ith);

    if (res_cst == 0) {
        return false;
    } else if (res_cst==2 && only_csp) {
        data_x->clear_fsbl_list(); // need to delete elements of fsbl_point_list since this branch is closed and they will not be needed later
        delete x_cell; // need to delete x_cell because not deleted in check cst.
        return false;
    } else if (data_x->pu != 1) {
        x_ctc[ith]->contract(x_cell->box);
        if(x_cell->box.is_empty()) {
            delete x_cell;
            return false;
        }
    }

    if(!only_csp) {
        //************* point evaluation ****************
        for (int i=0 ; i < nb_point ; i++) {
            Cell *x_copy = new Cell(*x_cell); // copy of the Cell and the y_heap
            bool found_feas_pt = get_feasible_point(x_copy, ith);
            if (found_feas_pt) { // we found a feasible point
                lsolve[ith]->nb_iter = choose_nbiter(true, false, x_cell);
                lsolve[ith]->prec_y = prec_y;
                lsolve[ith]->list_elem_max = 0; // no limit on heap size
                lsolve[ith]->visit_all = false; // no need to visit all leaves in midpoint
                bool res1;

                if (min_goal) {

                    res1 = eval_goal(x_copy, loup, ith);
                } else {
                    lsolve[ith]->nb_iter = choose_nbiter(true, false, x_cell);
                    lsolve[ith]->prec_y = prec_y;
                    lsolve[ith]->list_elem_max = 0; // no limit on heap size
                    lsolve[ith]->visit_all = false; // no need to visit all leaves in midpoint
                    //cout<<"a4"<<endl;
                    res1 = lsolve[ith]->optimize(x_copy, loup); // eval maxf(midp,heap_copy), go to minimum prec on y to get a thin enclosure
                    //cout<<"a5"<<endl;
                    DataMinMaxOpti * data_x_copy = &(x_copy->get<DataMinMaxOpti>());
                }

                lsolve[ith]->visit_all = visit_all; // reset visit all to initial value
                if (!res1) {
                    delete x_copy;
                } else {
                    //cout<<"test4"<<endl;
                    Interval ev = x_copy->get<DataMinMaxOpti>().fmax;
                    IntervalVector ysol = x_copy->get<DataMinMaxOpti>().y_heap->top()->box;
                    Vector sol = lsolve[ith]->best_point_eval.mid();
                    double new_loup = x_copy->get<DataMinMaxOpti>().fmax.ub();
                    //cout << "new loup : " << new_loup << " ; loup : " << loup << endl;
                    if (new_loup < loup) { // update best current solution
                        omp_set_lock(&bufferlock);
                        loup = new_loup; // TODO omp : mutex on loup
                        loup_changed = true; // TODO omp : mutex on loup_changed
                        loup_point = (x_copy->box.mid());
                        omp_unset_lock(&bufferlock);
                        if (trace) cout << "[mid]"  << " loup update " << loup  << " loup point  " << loup_point << endl;
                    }
                    
                    delete x_copy; // delete copy of the heap, no more use and it was not delete in light optim since res1 = 1
                }
            }
        }

        //************ evaluation of f(x,y_heap) *****************
        lsolve[ith]->prec_y = compute_min_prec(x_cell->box,false);
        lsolve[ith]->nb_iter = choose_nbiter(false, false, x_cell);
        lsolve[ith]->list_elem_max = compute_heap_max_size(data_x->y_heap->size(), false);

        bool res;
        
        if(min_goal) { // not a minmax problem, only a min problem-> suffices
            res = eval_goal(x_cell, loup, ith);
        } else {
            res = lsolve[ith]->optimize(x_cell, loup);
            double how_far = (data_x->fmax.lb()-uplo)/uplo;
            how_far = 2;
        }

        if(!res) {
            delete x_cell;
            return false;
        }
    }

    //***** if x_cell is too small ******************
    if(x_cell->box.max_diam() < prec) {
        //cout << "Min prec, box: " << x_cell->box << endl;
        int ctr_ok = 2;
        tmpbox = x_cell->box;
        ctr_ok = check_constraints(x_cell, true, ith);

        if (ctr_ok !=0 && !only_csp) {
            lsolve[ith]->nb_iter = choose_nbiter(true, false, x_cell);   // need to be great enough so the minimum precision on y is reached
            lsolve[ith]->prec_y = prec_y;
            lsolve[ith]->list_elem_max = 0; // no limit on heap size
            bool res;
            if(min_goal) {
                res = eval_goal(x_cell, loup, ith);
            } else {
                res = lsolve[ith]->optimize(x_cell, loup); // eval maxf(midp,heap_copy), go to minimum prec on y to get a thin enclosure
            }
            if(!res) {
                data_x->clear_fsbl_list(); // need to delete elements of fsbl_point_list since this branch is closed and they will not be needed later
                delete x_cell;
            } else {
                omp_set_lock(&bufferlock);
                update_uplo_of_epsboxes(data_x->fmax.lb());
                omp_unset_lock(&bufferlock);
                data_x->clear_fsbl_list(); // need to delete elements of fsbl_point_list since this branch is closed and they will not be needed later
                delete x_cell;
            }
        }
        return false;
    }

    // update optim data of the cell
    omp_set_lock(&bufferlock);
    buffer->cost1().set_optim_data(*x_cell, *x_sys[ith]);
    buffer->cost2().set_optim_data(*x_cell, *x_sys[ith]);
    buffer->push(x_cell);
    nb_cells++;
    omp_unset_lock(&bufferlock);
    return true;
}


bool OptimMinMax::eval_goal(Cell* x_cell, double loup, int ith) {

    (x_sys[ith]->goal)->backward(Interval(NEG_INFINITY, loup), x_cell->box);
    if(x_cell->box.is_empty()) { return false; }
    DataMinMaxOpti * data_x = &(x_cell->get<DataMinMaxOpti>());
    data_x->fmax = (x_sys[ith]->goal)->eval(x_cell->box);
    if(data_x->fmax.lb()>loup) { return false; }
    return true;
}


double OptimMinMax::compute_min_prec( const IntervalVector& x_box, bool csp) {

    if (!csp) {
        if(min_prec_coef == 0)
            return prec_y;

        double new_prec = x_box.max_diam()*min_prec_coef;
        return new_prec>prec_y?new_prec:prec_y;
    } else {
        if(min_prec_coef_csp == 0)
            return prec_y;

        double new_prec = x_box.max_diam()*min_prec_coef;
        return new_prec>prec_fa_y?new_prec:prec_fa_y;
    }
}


int OptimMinMax::choose_nbiter(bool midpoint_eval, bool csp, Cell* x_cell) {

    if (!midpoint_eval) {
        if(!csp) {
            DataMinMaxOpti * data_x = &(x_cell->get<DataMinMaxOpti>());
            data_x->nb_bisect+=iter;
            if(propag)
                return iter;
            else
                return data_x->nb_bisect;
        } else {
            DataMinMaxCsp * data_x = &(x_cell->get<DataMinMaxCsp>());
            data_x->nb_bisect+=iter_csp;
            if(propag)
                return iter_csp;
            else
                return data_x->nb_bisect;
        }
    }
    else { return 0; }// go to min prec
}


int OptimMinMax::compute_heap_max_size(int y_heap_size, bool csp) {
    if (!csp) {
        if(list_rate == 0 || y_heap_size>= list_elem_absolute_max) { // no constraint on list growth rate or max size already reached
            return list_elem_absolute_max;
        }
        return  y_heap_size+list_rate;
    } else {
        if(list_rate_csp == 0 || y_heap_size>= list_elem_absolute_max_csp) { // no constraint on list growth rate or max size already reached
            return list_elem_absolute_max_csp;
        }
        return  y_heap_size+list_rate_csp;
    }
}

bool OptimMinMax::get_feasible_point(Cell * elem, int ith) {
    elem->box = elem->box.mid(); // elem->box = elem->box.mid(); //get the box (x,mid(y))
    int res = check_constraints(elem, true, ith);  // cout<<"done, res = "<<res<<endl;
    if(res == 2) { return true; }
    return false;
}


int OptimMinMax::check_regular_ctr(const IntervalVector& box, int ith) {
    int res =2;
    for(int i=0 ; i < x_sys[ith]->nb_ctr ; i++) {
        Interval int_res = (x_sys[ith]->ctrs[i]).f.eval(box);
        if(int_res.lb() >= 0) {
            return 0;
        } else if(int_res.ub() >= 0) {
            res = 1;
        }
    }
    return res;
}


int OptimMinMax::check_fa_ctr(Cell* x_cell, bool midp, int ith) {

    DataMinMax * data_csp = &(x_cell->get<DataMinMaxCsp>());

    int res = 2;
    
    //cout << "data_csp = " << data_csp->pu << endl;
    if(data_csp->pu != 1) {

        fa_lsolve[ith]->nb_iter = choose_nbiter(midp, true, x_cell);
        //cout << "nb_iter = " << fa_lsolve[ith]->nb_iter << endl;
        if(midp) {
            fa_lsolve[ith]->visit_all = false; // no need to visit the heap in midp
            fa_lsolve[ith]->list_elem_max = 0;
            fa_lsolve[ith]->prec_y = prec_fa_y;
        }
        else {
            fa_lsolve[ith]->list_elem_max = compute_heap_max_size(data_csp->y_heap->size(), true);
            fa_lsolve[ith]->prec_y = compute_min_prec(x_cell->box, true);
        }
        //cout << "start fa_lsolve" << endl;
        bool ok = fa_lsolve[ith]->optimize(x_cell, 0);
        //cout << "end fa_lsolve" << endl;
        fa_lsolve[ith]->visit_all = visit_all_csp; // reset visit all to initial value
        if(!ok) {
            res = 0;
        } else {
            if(data_csp->y_heap->empty()) {
                data_csp->pu = 1;
                res = 2;
            } else if(data_csp->y_heap->top1()->get<OptimData>().pf.ub() < 0) {
                data_csp->pu = 1;
                data_csp->y_heap->flush();
                res = 2;
            }

            else res = 1;
        }
    } else {
        res = 2;
    }
    //if (!midp) cout << res << endl;
    return res;
}


int OptimMinMax::check_constraints(Cell * x_cell, bool midp, int ith) {

    int res_rctr  = 2;
    int res_factr = 2;

    DataMinMaxOpti * data_opt = &(x_cell->get<DataMinMaxOpti>());
    if(data_opt->pu != 1) {
        res_rctr = check_regular_ctr(x_cell->box, ith);
    }
    if(res_rctr == 2) {
        data_opt->pu = 1;
    } else if(res_rctr == 0) {
        if(!midp) { data_opt->clear_fsbl_list(); }
        delete x_cell;
        return 0;
    }
    if(fa_y_cst) { res_factr = check_fa_ctr(x_cell, midp, ith); }

    if(res_factr==0) { delete x_cell; }
    if(res_rctr == 2 &&  res_factr == 2) { // all ctr satisfied
        return 2;
    } else if(res_rctr == 0 || res_factr == 0) {
        return 0;
    } else {
        return 1;
    }
}


void OptimMinMax::show_yheap(Cell * x_cell, int ith) {

    DataMinMaxCsp * data_opt = &(x_cell->get<DataMinMaxCsp>());
    std::vector<Cell *> y_cells;
    Cell * y_cell;
    OptimData  *data_y;

    while(!data_opt->y_heap->empty()) {
        y_cell = data_opt->y_heap->pop1();
        data_y = &(y_cell->get<OptimData>());
        IntervalVector box(x_cell->box.size()+y_cell->box.size());
        box.put(0, x_cell->box);
        box.put(x_cell->box.size(), y_cell->box);
        y_cells.push_back(y_cell);
    }

    while(!y_cells.empty()) {
        y_cell = y_cells.back();
        fa_lsolve[ith]->check_already_in(y_cell, data_opt->y_heap);
        data_opt->y_heap->push(y_cell);
        y_cells.pop_back();
    }
}

}


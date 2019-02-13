//============================================================================
//                                  I B E X
// File        : ibex_LightOptimMinMax.cpp
// Author      : Dominique Monnet, Jordan Ninin
// License     : See the LICENSE file
// Created     : Oct 1, 2016
//============================================================================


#include "ibex_LightOptimMinMax.h"
#include "ibex_DataMinMax.h"
#include "ibex_NoBisectableVariableException.h"
#include "ibex_Timer.h"


using namespace std;

namespace ibex{


const double LightOptimMinMax::default_timeout = 2;
const double LightOptimMinMax::default_goal_abs_prec = 0;


LightOptimMinMax::LightOptimMinMax(NormalizedSystem& y_sys     , 
                                   Ctc&              ctc_xy    ,
                                   bool              csp_actif):
    timeout(default_timeout),
    ctc_xy(ctc_xy),
    xy_sys(y_sys),
    bsc(new LargestFirst()), prec_y(0), found_point(false), time(0),
    list_elem_max(0),
    nb_iter(0),
    csp_actif(csp_actif),
    goal_abs_prec(default_goal_abs_prec),
    best_point_eval(y_sys.box){};


LightOptimMinMax::~LightOptimMinMax() { delete bsc; }


void LightOptimMinMax::add_backtrackable(Cell& root, const IntervalVector& y_init,int critpr) {
    if (csp_actif) {root.add<DataMinMaxCsp>() ; }
    else           {root.add<DataMinMaxOpti>(); }

    Cell * y_cell = new Cell(y_init);
    y_cell->add<OptimData>();
    bsc->add_backtrackable(*y_cell);
    DataMinMax *data_x;

    if (csp_actif) { data_x = &(root.get<DataMinMaxCsp>()) ; }
    else           { data_x = &(root.get<DataMinMaxOpti>()); }

    data_x->y_heap_costf1.add_backtrackable(*y_cell);
    data_x->y_heap_costf2.add_backtrackable(*y_cell);
    data_x->y_heap->push(y_cell);
    data_x->y_heap->critpr = critpr;
}

bool LightOptimMinMax::optimize(Cell* x_cell, double loup) {

    delete_save_heap();

    found_point = false;
    DataMinMax *data_x;

    if (csp_actif ) { data_x = &(x_cell->get<DataMinMaxCsp>()) ; }
    else            { data_x = &(x_cell->get<DataMinMaxOpti>()); }

    DoubleHeap<Cell> *y_heap = data_x->y_heap;

    time = Timer::get_time();

    // ********** contract x_box with ctc_xy***************
    IntervalVector xy_box = xy_box_hull(x_cell->box);
    ctc_xy.contract(xy_box);

    if(xy_box.is_empty()) {
        return false;
    } else {
        x_cell->box &= xy_box.subvector(0,x_cell->box.size()-1);
    }

    std::vector<double> ub,lb,nbel,nbel_save;

    save_heap_ub = NEG_INFINITY;

    // *********** loop ********************
    try {
        int current_iter = 1;
        //cout << "start light loop " << endl;
        while(!stop_crit_reached(current_iter, y_heap,data_x->fmax)) {
            found_point  = false; 

            Cell * y_cell = y_heap->pop(); // we extract an element with critprob probability to take it according to the first crit
            current_iter++;
            //if (csp_actif) cout << "fa_lsolve3 : " << current_iter << endl;
            
            if((list_elem_max != 0 && ((y_heap->size() + heap_save.size())>list_elem_max)) || (y_cell->box.size())<prec_y) {
                bool res = handle_cell( x_cell, y_cell,loup);
                if (!res) { return false; }// x_cell has been deleted
            } else {
                try {
                    std::pair<Cell*,Cell*> subcells_pair=bsc->bisect_cell(*y_cell);// bisect tmp_cell into 2 subcells
                    //if (csp_actif) cout << "y_cell = " << *y_cell << endl;
                    delete y_cell;
                    bool res = handle_cell( x_cell, subcells_pair.first,loup);
                    if (!res) { // x_cell has been deleted
                        delete subcells_pair.second;
                        return false;
                    }
                    res = handle_cell( x_cell, subcells_pair.second,loup);
                    if (!res) {return false; } // x_cell has been deleted
                }
                catch (NoBisectableVariableException& ) {
                    bool res = handle_cell(x_cell,y_cell,loup);
                    if(!res) { return false; }
                }
            }
            Timer::check(time+timeout);
        }
    }
    catch (TimeOutException& ) { }

    if (visit_all == true) {
        while(!y_heap->empty()) {
            Cell * y_cell = y_heap->pop();
            bool res = handle_cell(x_cell,y_cell,loup,true);
            if (!res) { return false; }
        }
    }

    fill_y_heap(*y_heap);

    if(y_heap->empty()){
        if(csp_actif) { return true; }// sic case: empty set y means that the set defined by the constraints g(x,y)<0 is empty, and therefore that the sic is respected over the empty set 
        else { return false; }// minmax case: empty set means ????? suppose that x can be discarded????
    }

    double new_fmax_ub = y_heap->top1()->get<OptimData>().pf.ub(); // get the upper bound of max f(x,y_heap)
    double new_fmax_lb = y_heap->top2()->get<OptimData>().pf.lb(); // get the lower bound of max f(x,y_heap)

    if (new_fmax_ub < new_fmax_lb) {
        ibex_error("ibex_LightOptimMinMax: error, please report this bug.");
    }

    data_x->fmax &= Interval(new_fmax_lb, new_fmax_ub);

    if(  data_x->fmax.is_empty() || data_x->fmax.lb() > loup) {
        return false;
    }

    best_point_eval = xy_box.mid();
    for(int i=0;i<xy_sys.box.size()-x_cell->box.size();i++) {
        best_point_eval[x_cell->box.size()+i] = y_heap->top1()->box[i].mid();
    }
    return true;
}


bool LightOptimMinMax::stop_crit_reached(int current_iter,DoubleHeap<Cell> * y_heap,const Interval& fmax) {

    if(nb_iter !=0 && current_iter>=nb_iter) { return true; } // nb_iter ==  0 implies minimum precision required (may be mid point x case)
    if(y_heap->size() == 0) { return true; }
    if(csp_actif && (y_heap->top1()->get<OptimData>().pf.ub() < 0) ) { return true; }
    return false;
}


bool LightOptimMinMax::handle_cell( Cell* x_cell,Cell*  y_cell,double loup,bool no_stack) {

    IntervalVector xy_box =init_xy_box(x_cell->box,y_cell->box);
    // recuperer les data
    DataMinMax *data_x;

    if (csp_actif) { data_x = &(x_cell->get<DataMinMaxCsp>()); }
    else { data_x = &(x_cell->get<DataMinMaxOpti>()); }

    OptimData  *data_y = &(y_cell->get<OptimData>());

    if(data_y->pu!=1) { // Check constraints
        if(handle_constraint(data_y, xy_box, y_cell->box)) {
            delete y_cell;
            return true;
        }
    } else { handle_cstfree(xy_box,y_cell); }

    /********************************************************************************/
    IntervalVector mid_y_box = get_feasible_point(x_cell,y_cell);

    if (!(mid_y_box.is_empty())) {
        Interval midres = eval_all(xy_sys.goal,mid_y_box);
        if ( loup < midres.lb() ) {
            delete y_cell;
            return false; // no need to go further, x_box does not contains the solution
        } else if(midres.lb()>data_y->pf.lb()) {
            data_y->pf   &= Interval(midres.lb(),POS_INFINITY);
            if(data_x->best_sol!=NULL) { delete data_x->best_sol; }
            data_x->best_sol = new IntervalVector(mid_y_box);
            if (data_x->fmax.lb() < midres.lb() ) {
                found_point = true;
                data_x->fmax &= Interval(midres.lb(),POS_INFINITY);; // yes we found a feasible solution for all x
            }
        }
    }

    //************ part below add a contraction w.r.t f(x,y)<best_max, this part may not be efficient on every problem ******************************

    if(data_y->pu == 1) {
        IntervalVector xy_box_mem(xy_box);
        xy_sys.goal->backward(Interval(NEG_INFINITY,loup),mid_y_box);

        if(xy_box.is_empty()) {
            delete y_cell;
            return false;
        }
        for(int i=x_cell->box.size();i<xy_box.size();i++) { // y contracted => E y, for all x f(x,y)>loup, x deleted
            if(xy_box[i] != xy_box_mem[i]) {
                delete y_cell;
                return false;
            }
        }

        for (int k=0; k<x_cell->box.size(); k++) {
            x_cell->box[k] &= xy_box[k];
        }
    }
    //********************************************
    if(eval_all(xy_sys.goal,xy_box).ub() > data_y->pf.ub() + 10e-13) {
       cout<<" ************************** CRITICAL ISSUE *******************"<<endl;
       cout<<" get worst upper bound, should not happen due to monotonicity of ifunc"<<endl;
       cout << xy_box << endl;
       cout << eval_all(xy_sys.goal,xy_box).ub() - data_y->pf.ub() << endl;

       cout<<"***************************************************************"<<endl;
    }

    data_y->pf &= eval_all(xy_sys.goal,xy_box);
    if( data_y->pf.is_empty() || data_x->fmax.lb() > data_y->pf.ub()) {  // y_box cannot contains max f(x,y)
        delete y_cell;
        return true;
    }

    if((data_y->pf.lb() > loup) && (data_y->pu == 1)) {
        delete y_cell;
        return false; // no need to go further, x_box does not contains the solution
    }

    //*************************************************
    if (y_cell->box.max_diam()<prec_y) {
        save_heap_ub = save_heap_ub<data_y->pf.ub()?data_y->pf.ub():save_heap_ub;
        heap_save.push_back(y_cell);
    } else {
        if(!no_stack) {
            data_x->y_heap->push(y_cell);
        } else {
            heap_save.push_back(y_cell);
        }
    }
    return true;
}


bool LightOptimMinMax::handle_constraint(OptimData  *data_y, IntervalVector & xy_box,IntervalVector & y_box) {

    switch(check_constraints(xy_box)) {
        case 2:  { data_y->pu = 1; break; }
        case 0:  { return true; }
        default: { break; }
    }

    if(data_y->pu != 1)  {
        ctc_xy.contract(xy_box);
        if (xy_box.is_empty()) {
            return true;
        } else {
            for (int k=0; k<y_box.size(); k++) {
                y_box[k] &= xy_box[xy_box.size()-y_box.size()+k];
            }
        }
    }
    return false;
}


bool LightOptimMinMax::handle_cstfree(IntervalVector& xy_box, Cell * const y_cell) {

    IntervalVector grad(xy_box.size());
    xy_sys.goal->gradient(xy_box,grad);
    for(int i = xy_box.size() - y_cell->box.size() ; i < xy_box.size() ; i++) {
        if(grad[i].lb()>0) { (xy_box)[i] = Interval((xy_box)[i].ub()); }
        if(grad[i].ub()<0) { (xy_box)[i] = Interval((xy_box)[i].lb()); }
    }
    return true;
}


IntervalVector LightOptimMinMax::get_feasible_point(Cell * x_cell,Cell * const y_cell) {

    IntervalVector mid_y_box = get_mid_y(x_cell->box,y_cell->box); // get the box (x,mid(y))
    if((y_cell->get<OptimData>().pu != 1)) { // constraint on xy exist and is not proved to be satisfied
        int res = check_constraints(mid_y_box);
        if(res == 0 ||res == 1) { return IntervalVector(1,Interval::EMPTY_SET); }
    }
    return mid_y_box;
}


int LightOptimMinMax::check_constraints(const IntervalVector& xy_box) {
    int res = 2;

    for(int i = 0 ; i < xy_sys.nb_ctr ; i++) {
        Interval int_res =  eval_all(&(xy_sys.ctrs[i].f), xy_box);
        if(int_res.lb() > 0) { 
            return 0;
        } else if(int_res.ub() >= 0) {
            res = 1;
        }
    }
    return res;
}


IntervalVector LightOptimMinMax::get_mid_y(const IntervalVector& x_box, const IntervalVector& y_box) { // returns the cast of box x and mid of box y

    IntervalVector res(x_box.size() + y_box.size());
    for (int i = 0 ; i < x_box.size() ; i++) { res[i] = x_box[i];}
    for (int i = 0 ; i < y_box.size() ; i++) { res[x_box.size()+i] = y_box[i].mid(); }
    return res;
}


IntervalVector LightOptimMinMax::init_xy_box(const IntervalVector& x_box, const IntervalVector & y_box) {
    IntervalVector res(x_box.size() + y_box.size());
    for (int k = 0 ; k < x_box.size() ; k++) { res[k] = x_box[k]; }
    for (int k = 0 ; k < y_box.size() ; k++) { res[k+x_box.size()] = y_box[k]; }
    return res;
}


IntervalVector LightOptimMinMax::xy_box_hull(const IntervalVector& x_box) {
    IntervalVector res(xy_sys.nb_var);
    for(int k = 0 ; k<x_box.size() ; k++) { res[k] = x_box[k]; }
    for(int k = x_box.size();k < xy_sys.nb_var ; k++) // update current y box in the xy_box
        res[k] = xy_sys.box[k];
    return res;
}


void LightOptimMinMax::fill_y_heap(DoubleHeap<Cell>& y_heap) {
    Cell* tmp_cell;
    while(!heap_save.empty()) { // push all boxes of heap_save in y_heap
        tmp_cell = heap_save.back();
        y_heap.push(tmp_cell);
        heap_save.pop_back();
    }
}


void LightOptimMinMax::delete_save_heap() {
    while(!heap_save.empty()) {
        delete heap_save.back();
        heap_save.pop_back();
    }
}


Interval LightOptimMinMax::eval_all(Function* f,const IntervalVector& box) {
    Interval eval = f->eval(box);
    return eval;
}

bool LightOptimMinMax::check_already_in(Cell * const y_cell, DoubleHeap<Cell> * y_heap) {
    vector<Cell *> stack;
    while(!y_heap->empty()) {
        Cell * c = y_heap->pop1();
        stack.push_back(c);
        if(y_cell == c) {
            cout<<" ***************** CRITICAL ERROR****************"<<endl;
            cout<<"   insert same cell in the heap, pointer =  "<<y_cell<<",  box: "<<y_cell->box<<endl;
            cout<<" ************************************************"<<endl;
        }
    }
    while(!stack.empty()) {
        y_heap->push(stack.back());
        stack.pop_back();
    }

}

}



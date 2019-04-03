#include "ibex_SortedHeap.h"


Interval evaluate_surface(SortedHeap& w_heap)
{
    // Evaluate the surface

    double lb = 0;
    double ub = 0;
    IntervalVector res_box = IntervalVector(2);

    while(!w_heap.is_empty())
    {
        res_box = w_heap.pop();
        lb += (res_box[0]).lb() * (res_box[1]).diam();
        ub += (res_box[0]).ub() * (res_box[1]).diam();
    }

    cout << Interval(lb,ub) << endl;
    return Interval(lb,ub);
}


Interval eval_integral(const Function& f, const IntervalVector& box, int max_w_cell, Interval& w)
{
    int n = box.size();

    IntervalVector box_with_w = box;

    box_with_w.resize(n+1); 

    box_with_w.put(n,IntervalVector(1,w));

    SortedHeap w_heap = SortedHeap();
    
    Interval res = f.eval(box_with_w);
    IntervalVector first_res = IntervalVector(2);
    first_res[0] = res;
    first_res[1] = w;

    w_heap.push(first_res);
    
    for (int iter = 0 ; iter < max_w_cell ; iter++)
    {

        IntervalVector box_max = w_heap.pop();        
        std::pair<IntervalVector,IntervalVector> p = box_max.bisect(1,0.5);

        IntervalVector box_res_1_init = p.first;
        IntervalVector box_res_2_init = p.second;

        IntervalVector box_1 = box_with_w;
        IntervalVector box_2 = box_with_w;
        
        Interval w_1 = box_res_1_init[1];
        Interval w_2 = box_res_2_init[1];
        box_1.put(n,IntervalVector(1,w_1));
        box_2.put(n,IntervalVector(1,w_2));

        Interval inter_res_1 = f.eval(box_1);
        Interval inter_res_2 = f.eval(box_2);

        IntervalVector res_1 = IntervalVector(2);
        IntervalVector res_2 = IntervalVector(2);

        res_1[0] = inter_res_1; res_1[1] = w_1;
        res_2[0] = inter_res_2; res_2[1] = w_2;

        w_heap.push(res_1);
        w_heap.push(res_2);

        // SortedHeap new_heap = w_heap;
        // evaluate_surface(new_heap);
    } 

    return evaluate_surface(w_heap);
}



int main(int argc, char* argv[]) {

    int max_w_cell = 10000;
    Interval w_range = Interval(0,10);

    Variable x(4),w;
    IntervalVector x_ini(4);

    x_ini[0] = Interval(3,3);
    x_ini[1] = Interval(3,3);
    x_ini[2] = Interval(3,3); 
    x_ini[3] = Interval(3,3);

    Function f = Function(x,w, pow(x[0],2)*w*x[1]+ x[3]/(w+x[2]*x[1]));

    cout << f << endl;

    time_t start, end, diff;
    start = clock();
    Interval iv = eval_integral(f, x_ini, max_w_cell, w_range);
    end   = clock();
    cout << "difftime = " << ((double) (end - start)) / CLOCKS_PER_SEC << endl;

    cout << iv << endl;
    return 1;

}
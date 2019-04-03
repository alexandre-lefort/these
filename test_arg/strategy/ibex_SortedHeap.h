#ifndef __SORTED_HEAP__
#define __SORTED_HEAP__


#include "ibex_IntervalVector.h"
#include <vector>
#include <iostream>


using namespace std;
using namespace ibex;


class SortedHeap {

private:

    std::vector<IntervalVector> heap;

public:
    SortedHeap();
    
    void push(IntervalVector& new_box);
    IntervalVector pop();
    bool is_empty();
};


#endif
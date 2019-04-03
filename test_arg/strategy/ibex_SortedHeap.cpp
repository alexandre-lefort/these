#include "ibex_SortedHeap.h"


SortedHeap::SortedHeap() {};



bool compare(const IntervalVector& t1, const IntervalVector& t2)
{
    return (t1.volume() <= t2.volume());
}


void SortedHeap::push(IntervalVector& new_box)
{
    int n = heap.size();
    int aux = 0;

    if (n == 0) {
        heap.push_back(new_box);
        return;
    }
    
    vector<IntervalVector>::iterator itt = std::upper_bound(heap.begin(),heap.end(), new_box, compare);

    aux = distance(heap.begin(),itt);

    if (aux == (n-1)) {
        heap.push_back(new_box);
        return;
    }

    heap.insert(heap.begin()+aux,new_box);

}


IntervalVector SortedHeap::pop()
{
    IntervalVector res = heap.back();
    heap.pop_back();
    return res;
};


bool SortedHeap::is_empty()
{
    return (heap.size() < 1);
};
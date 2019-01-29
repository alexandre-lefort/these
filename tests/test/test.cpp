#include "omp.h"
#include <cstdio>
#include <iostream>

int main()
{
omp_set_num_threads(8);
#pragma omp parallel
{
    std::cout << omp_get_num_threads() << std::endl;
    int ID = omp_get_thread_num();
    printf(" hello(%d)", ID);
    printf(" world(%d) \n", ID);
}
return 1;
}

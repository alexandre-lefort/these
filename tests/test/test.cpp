#include "omp.h"
#include <cstdio>
#include <iostream>

int main()
{
omp_set_num_threads(16);
std::cout << "number of threads :" <<  omp_get_num_threads() << std::endl;
#pragma omp parallel
{   
	std::cout << "number of threads :" <<  omp_get_num_threads() << std::endl;
    
    printf(" hello(%d)", ID);
    printf(" world(%d) \n", ID);
}
std::cout << "number of threads :" <<  omp_get_num_threads() << std::endl;
return 1;
}

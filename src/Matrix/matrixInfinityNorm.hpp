#ifndef matrixInfinityNorm_HPP
#define matrixInfinityNorm_HPP
#include "Matrix.hpp"

#define MAX(a, b)(a>b?a:b)
#define ABS(a)(a>0?a:(-a))

float matrixInfinityNorm(Matrix &A)
{
    float norm=0.0,sum;
    unsigned i,j;

    for(i=0;i<A.NR;i++)
    {
        sum =   0.0;
        for(j=0;j<A.NC;j++)
            sum +=ABS(A[i][j]);
        norm    =   MAX(norm,sum);
    }

    return norm;
}



#endif

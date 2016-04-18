#include "src/includes.hpp"

#define NR 3
#define NC 3

int main()
{
    Matrix A(NR,NC);
    unsigned i,j;

    for(i=0;i<NR;i++)
        for(j=0;j<NC;j++)
            A[i][j] =   -(1.0*i*i+1.0*j*j);

    matrixDisplay(A);
    printf("Norm of the Matrix=%6.2f\n",matrixInfinityNorm(A));
    return 0;
}

#include "src/includes.hpp"

#define N 1000

int main()
{
    Matrix A(N,N);
    Matrix R(N,N);
    unsigned i,j;

    for(i=0;i<N;i++)
        for(j=0;j<N;j++)
        {
            A[i][j] =   0.01;
        }

    matrixExponential5_serial(A,R);
    //printf("**A**\n");
    //matrixDisplay(A);
    printf("*****\n");
    printf("**R**\n");
    matrixDisplay(R);
    printf("*****\n");
    return 0;
}

#include "src/includes.hpp"

#define NR 3
#define NC 3

int main()
{
    Matrix A(NR,NC);
    Matrix B(NR,3);
    Matrix X(NR,3);
    unsigned i,j;

    for(i=0;i<NC;i++)
        for(j=0;j<NC;j++)
        {
            A[i][j] =   i+j+2;
            B[i][j] =   1;
        }

    A[0][2] =   63;
    A[2][2] =   163;
    B[0][0] =   3;
    B[1][0] =   4;
    B[2][0] =   1;

    solveAXB(A,X,B);

    matrixDisplay(X);
    printf("*****\n");
    solveAXB_openmp(A,X,B);

    matrixDisplay(X);
    printf("*****\n");
    return 0;
}

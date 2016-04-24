#include "src/includes.hpp"

#define NR 3
#define NC 3

int main()
{
    Matrix A(NR,NC);
    Matrix B(NR,NC);
    unsigned i,j;

    for(i=0;i<NC;i++)
        for(j=0;j<NC;j++)
            A[i][j] =   i*i +j*j*j*j*j*j;

    A[0][2] =   63;

    printf("**A**\n");
    matrixDisplay(A);
    printf("*****\n");
    matrixInverse(A,B);
    printf("**B**\n");
    matrixDisplay(B);
    printf("*****\n");
    return 0;
}

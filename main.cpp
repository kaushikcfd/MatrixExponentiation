#include "src/includes.hpp"

<<<<<<< HEAD
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
=======
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

    printf("**A**\n");
    matrixDisplay(A);
    printf("*****\n");
    printf("**B**\n");
    matrixDisplay(B);
    printf("*****\n");
    printf("**X**\n");
    matrixDisplay(X);
>>>>>>> origin/master
    printf("*****\n");
    return 0;
}

#include "Matrix.hpp"
#include "matrixInfinityNorm.hpp"
#include "matrixMatrixMultiply.hpp"
#include "solveAXB.hpp"
#include "matrixDisplay.hpp"
#include <cblas.h>
#include <cstring>
#include <cmath>/*This would be used for calculating the log to the base 2 for finding the squaring coefficient.*/

#ifndef matrixExponential_serial_HPP
#define matrixExponential_serial_HPP
#define MAX(a, b)(a>b?a:b)

void matrixExponential_serial(Matrix &A, Matrix &Exp)
{
    int sigma,i,j;
    float norm,coeff;/*`sigma` is for the scaling down purposes. `norm` is the infinty norm which would be helpful for calculating the scaling down coefficient. And `coeff` is $2^sigma$*/
    if((A.NC!=A.NR)||(Exp.NC!=Exp.NR)||(Exp.NC!=Exp.NR))/*This condition is checking whether the matrices given are compatible to get the Matrix Exponential. The checks performed are 1.) Both are square matrices and 2.) Of the same order. */
    {
        fprintf(stderr,"The given Matrices were not of square type hence it was not possible to compute the Maatrix Exponential. EXITING!.\n");
        return ;
    }
    unsigned N =    A.NC;/*Getting the order of the matrix.*/
    Matrix R(N,N),R_2(N,N);/*This is the Pade Approximant.*/
    Matrix B(N,N),B_2(N,N),B_3(N,N);
    Matrix Identity(N,N);
    Matrix Num(N,N);/*This would form the numerator of the Pade Approximant.*/
    Matrix Den(N,N);/*This would form the denominator of the Pade Approximant.*/
    memcpy(&B[0][0],&A[0][0],N*N*sizeof(float));/*Copying it into a dummy variable so that the structure of given matrix does not change.*/
    /*Finding the scaling factor.*/
    norm    =   matrixInfinityNorm(A);
    sigma   =   MAX(0,(int)log2(norm));
    coeff   =   pow(2,-sigma);
    cblas_sscal((N*N),coeff,&B[0][0],1);
    matrixMatrixMultiply_fast(B,B,B_2);
    matrixMatrixMultiply_fast(B_2,B,B_3);
    /*Initializing the Identity Matrix, Numerator Matrix and Denominator matrices.*/
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            Num[i][j]       =   0.0;
            Den[i][j]       =   0.0;
            Identity[i][j]  =   0.0;
        }
        Identity[i][i]      =   1.0;
    }
    cblas_saxpy((N*N),120.0,&Identity[0][0],1,&Num[0][0],1);
    cblas_saxpy((N*N),60.0,&B[0][0],1,&Num[0][0],1);
    cblas_saxpy((N*N),12.0,&B_2[0][0],1,&Num[0][0],1);
    cblas_saxpy((N*N),1.0,&B_3[0][0],1,&Num[0][0],1);

    cblas_saxpy((N*N),120.0,&Identity[0][0],1,&Den[0][0],1);
    cblas_saxpy((N*N),-60.0,&B[0][0],1,&Den[0][0],1);
    cblas_saxpy((N*N),12.0,&B_2[0][0],1,&Den[0][0],1);
    cblas_saxpy((N*N),-1.0,&B_3[0][0],1,&Den[0][0],1);

    solveAXB(Den,R,Num);
    for(i=0;i<sigma;i++)
    {
        matrixMatrixMultiply(R,R,R_2);
        memcpy(&R[0][0],&R_2[0][0],N*N*sizeof(float));
    }
    //TODO: Check this once.
    memcpy(&Exp[0][0],&R[0][0],N*N*sizeof(float));
    return ;
}

void matrixExponential5_serial(Matrix &A, Matrix &Exp)
{
    int sigma,i,j;
    float norm,coeff;/*`sigma` is for the scaling down purposes. `norm` is the infinty norm which would be helpful for calculating the scaling down coefficient. And `coeff` is $2^sigma$*/
    if((A.NC!=A.NR)||(Exp.NC!=Exp.NR)||(Exp.NC!=Exp.NR))/*This condition is checking whether the matrices given are compatible to get the Matrix Exponential. The checks performed are 1.) Both are square matrices and 2.) Of the same order. */
    {
        fprintf(stderr,"The given Matrices were not of square type hence it was not possible to compute the Maatrix Exponential. EXITING!.\n");
        return ;
    }
    unsigned N =    A.NC;/*Getting the order of the matrix.*/
    Matrix R(N,N),R_2(N,N);/*This is the Pade Approximant.*/
    Matrix B(N,N),B_2(N,N),B_3(N,N),B_4(N,N),B_5(N,N);
    Matrix Identity(N,N);
    Matrix Num(N,N);/*This would form the numerator of the Pade Approximant.*/
    Matrix Den(N,N);/*This would form the denominator of the Pade Approximant.*/
    memcpy(&B[0][0],&A[0][0],N*N*sizeof(float));/*Copying it into a dummy variable so that the structure of given matrix does not change.*/
    /*Finding the scaling factor.*/
    norm    =   matrixInfinityNorm(A);
    sigma   =   MAX(0,(int)log2(norm));
    coeff   =   pow(2,-sigma);
    cblas_sscal((N*N),coeff,&B[0][0],1);
    matrixMatrixMultiply_fast(B,B,B_2);
    matrixMatrixMultiply_fast(B_2,B,B_3);
    matrixMatrixMultiply_fast(B_3,B,B_4);
    matrixMatrixMultiply_fast(B_4,B,B_5);

    /*Initializing the Identity Matrix, Numerator Matrix and Denominator matrices.*/
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            Num[i][j]       =   0.0;
            Den[i][j]       =   0.0;
            Identity[i][j]  =   0.0;
        }
        Identity[i][i]      =   1.0;
    }
    cblas_saxpy((N*N),30240.0,&Identity[0][0],1,&Num[0][0],1);
    cblas_saxpy((N*N),15120.0,&B[0][0],1,&Num[0][0],1);
    cblas_saxpy((N*N),3360.0,&B_2[0][0],1,&Num[0][0],1);
    cblas_saxpy((N*N),420.0,&B_3[0][0],1,&Num[0][0],1);
    cblas_saxpy((N*N),30.0,&B_4[0][0],1,&Num[0][0],1);
    cblas_saxpy((N*N),1.0,&B_5[0][0],1,&Num[0][0],1);

    cblas_saxpy((N*N),30240.0,&Identity[0][0],1,&Den[0][0],1);
    cblas_saxpy((N*N),-15120.0,&B[0][0],1,&Den[0][0],1);
    cblas_saxpy((N*N),3360.0,&B_2[0][0],1,&Den[0][0],1);
    cblas_saxpy((N*N),-420.0,&B_3[0][0],1,&Den[0][0],1);
    cblas_saxpy((N*N),30.0,&B_4[0][0],1,&Den[0][0],1);
    cblas_saxpy((N*N),-1.0,&B_5[0][0],1,&Den[0][0],1);

    solveAXB(Den,R,Num);
    for(i=0;i<sigma;i++)
    {
        matrixMatrixMultiply_fast(R,R,R_2);
        memcpy(&R[0][0],&R_2[0][0],N*N*sizeof(float));
    }
    //TODO: Check this once.
    memcpy(&Exp[0][0],&R[0][0],N*N*sizeof(float));
    return ;
}

#endif

#include <lapacke.h>
#include <cstring>
#include <cstdio>
#include "Matrix.hpp"
#include "omp.h"

#ifndef SolveAXB_HPP
#define SolveAXB_HPP

#define MOD(a)(a>0?a:(-a))

void solveAXB(Matrix &A, Matrix &X, Matrix &B)
{
    if(((A.NR)!=(A.NC))||((X.NC)!=(B.NC))||((X.NR)!=(B.NR))||((X.NR)!=(A.NC)))
    {
        printf("The sizes of the given matrices are not right to perform Gaussian Elimination.\n");
        return ;
    }

    unsigned N      =   A.NC;
    unsigned NRHS   =   B.NC;

    float *a    = new float[N*N];
    float *b    = new float[N*NRHS];

    memcpy(a,&A[0][0],N*N*sizeof(float));/*Copying the matrix A into a dummy variable, so that the input matrix does not undergo ay change.*/
    memcpy(b,&B[0][0],N*NRHS*sizeof(float));/*Copying the matrix B into a dummy variable, so that the input matrix does not undergo ay change.*/
    int ipiv[N];
    int info;

    info    =   LAPACKE_sgesv(LAPACK_ROW_MAJOR,N,NRHS,a,N,ipiv,b,NRHS);/*Implementing the LAPACK subroutine for finding the solution.*/
    if(info!=0)
        fprintf(stderr,"The Linear solve `AX=B` was unsuccesful. Please Check the inputs\n");

    memcpy(&X[0][0],b,N*NRHS*sizeof(float));

    delete[] a;
    delete[] b;
    return ;
}

void solveAXB_openmp(Matrix &A, Matrix &X, Matrix &B)
{
    if(((A.NR)!=(A.NC))||((X.NC)!=(B.NC))||((X.NR)!=(B.NR))||((X.NR)!=(A.NC)))
    {
        printf("The sizes of the given matrices are not right to perform Gaussian Elimination.\n");
        return ;
    }

    unsigned N      =   A.NC;
    unsigned NRHS   =   B.NC;
    unsigned i;

    float *a    = new float[N*N];
    float *b    = new float[N*NRHS];

    memcpy(a,&A[0][0],N*N*sizeof(float));/*Copying the matrix A into a dummy variable, so that the input matrix does not undergo ay change.*/
    memcpy(b,&B[0][0],N*NRHS*sizeof(float));/*Copying the matrix B into a dummy variable, so that the input matrix does not undergo ay change.*/
    int ipiv[N];
    int info=0;

    # pragma omp parallel for private(i) shared (A,a,b)
    for(i=0;i<NRHS;i++)
    {
        info    =   LAPACKE_sgesv(LAPACK_ROW_MAJOR,N,1,a,N,ipiv,&b[i],NRHS);/*Implementing the LAPACK subroutine for finding the solution.*/
        memcpy(a,&A[0][0],N*N*sizeof(float));
        info+=MOD(info);
    }

    if(info!=0)
        fprintf(stderr,"The Linear solve `AX=B` was unsuccesful. Please Check the inputs\n");

    memcpy(&X[0][0],b,N*NRHS*sizeof(float));

    delete[] a;
    delete[] b;
    return ;
}



#endif

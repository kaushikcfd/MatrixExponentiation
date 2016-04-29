#ifndef matrixMatrixMultiply_HPP
#define matrixMatrixMultiply_HPP
#include "Matrix.hpp"
#include <lapacke.h>
#include <cblas.h>
#include <cstdio>
#include "omp.h"

void matrixMatrixMultiply( Matrix &A,  Matrix &B, Matrix &C )
{
    unsigned p  =   A.NR;
    unsigned q  =   A.NC  ;
    unsigned r  =   B.NR;
    unsigned s  =   B.NC;
    if((r!=q)||(C.NR!=p)||(C.NC!=s))
    {
        printf("The matrices do not commute. Exiting!\n");
        return ;
    }
    unsigned i,j,k;/*These variables act as counters for the loop.*/

    for(i=0;i<p;i++)
	    for(j=0;j<s;j++)
            {
                C[i][j] =   0.0;
                for(k=0;k<q;k++)
                    C[i][j] += (A[i][k]*B[j][k]);
            }

    return ;
}

void matrixMatrixMultiply_fast( Matrix &A,  Matrix &B, Matrix &C )
{
    unsigned p  =   A.NR;
    unsigned q  =   A.NC  ;
    unsigned r  =   B.NR;
    unsigned s  =   B.NC;
    if((r!=q)||(C.NR!=p)||(C.NC!=s))
    {
        printf("The matrices do not commute. Exiting!\n");
        return ;
    }
    cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,p,s,r,1.0,A[0],q,B[0],s,0.0,C[0],s);

    return ;
}

void matrixMatrixMultiply_openmp( Matrix &A,  Matrix &B, Matrix &C )
{
    unsigned p  =   A.NR;
    unsigned q  =   A.NC  ;
    unsigned r  =   B.NR;
    unsigned s  =   B.NC;
    unsigned i;
    if((r!=q)||(C.NR!=p)||(C.NC!=s))
    {
        printf("The matrices do not commute. Exiting!\n");
        return ;
    }

    # pragma omp parallel for private(i) shared (A,B,C)
    for(i=0;i<s;i++)
        cblas_sgemv(CblasRowMajor,CblasNoTrans,p,q,1.0,&A[0][0],q,&B[0][i],s,0.0,&C[0][i],s);

    return ;
}

#endif

#include <lapacke.h>
#include <cstring>
#include <cstdio>
#ifndef matrixInverse_HPP
#define matrixInverse_HPP

#include "Matrix.hpp"

void matrixInverse(Matrix &A, Matrix &Ainv)
{
    if((A.NC!=A.NR)||(Ainv.NC!=Ainv.NR)||(A.NC!=Ainv.NR))/*This condition is checking whether the matrices given are compatible to get the inverse. The checks performed are 1.) Both are square matrices and 2.) Of the same order. */
    {
        fprintf(stderr,"The given Matrices were not of square type hence it was not possible to compute the inverse. EXITING!.\n");
        return ;
    }
    unsigned N  =   A.NC;/*Setting `N` as the order of the matrix.*/

    memcpy(&Ainv[0][0],&A[0][0],N*N*sizeof(float));/*Making it ready for using **BLAS sgetri**.*/
    int info;
    int ipiv[N];

    info =  LAPACKE_sgetrf(LAPACK_ROW_MAJOR,N,N,&Ainv[0][0],N,ipiv);/*Getting the pivotization array ipiv[] for the sake of `sgetri`.*/
    if(info==0)
        info =  LAPACKE_sgetri(LAPACK_ROW_MAJOR,N,&Ainv[0][0],N,ipiv);/*Using the inverse finding function of the LAPACKE package.*/
    if(info!=0)
        fprintf(stderr,"The inverse of the matrix was unsuccesful.\n");/*Thrwoing an error when the inverse of the given matrix was not possible.*/
    return ;
}

#endif
/*END OF FILE*/

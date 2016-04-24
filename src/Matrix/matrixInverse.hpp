#include <lapacke.h>
#include <cstring>
#include <cstdio>
#ifndef matrixInverse_HPP
#define matrixInverse_HPP

#include "Matrix.hpp"

void matrixInverse(Matrix &A, Matrix &Ainv)
{
    if((A.NC!=A.NR)||(Ainv.NC!=Ainv.NR)||(A.NC!=Ainv.NR))
    {
        fprintf(stderr,"The given Matrices were not of square type hence it was not possible to compute the inverse. EXITING!.\n");
        return ;
    }
    unsigned N  =   A.NC;

    memcpy(&Ainv[0][0],&A[0][0],N*N*sizeof(float));
    int info;
    int ipiv[N];

    info =  LAPACKE_sgetrf(LAPACK_ROW_MAJOR,N,N,&Ainv[0][0],N,ipiv);
    if(info==0)
        info =  LAPACKE_sgetri(LAPACK_ROW_MAJOR,N,&Ainv[0][0],N,ipiv);
    if(info!=0)
        fprintf(stderr,"The inverse of the matrix was unsuccesful.\n");
    return ;
}

#endif

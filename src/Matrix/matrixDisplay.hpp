#ifndef matrixDisplay_HPP
#define matrixDisplay_HPP
#include "Matrix.hpp"
#include <cstdio>
void matrixDisplay(Matrix &A)
{
    unsigned i,j;
    for(i=0;i<A.NR;i++)
    {
        for(j=0;j<A.NC;j++)
            printf("%6.2f\t",A[i][j]);
        printf("\n");
    }
    return ;
}

void matrixDisplay4(Matrix &A)
{
    unsigned i,j;
    for(i=0;i<A.NR;i++)
    {
	for(j=0;j<A.NC;j++)
	    printf("%6.4f\t",A[i][j]);
	printf("\n");
    }
    return ;
}

#endif

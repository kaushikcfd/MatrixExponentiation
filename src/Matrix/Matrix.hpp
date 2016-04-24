#include <cstdio>
#ifndef Matrix_HPP
#define Matrix_HPP
class Matrix
{
public:
    float *Values;
    unsigned NR;
    unsigned NC;

    Matrix(unsigned , unsigned );
    ~Matrix();
    float* operator[](unsigned );
    void display();
};

Matrix::Matrix( unsigned m, unsigned n)
{
    NR  =   m;
    NC  =   n;
    Values  =   new float[NR*NC];
}

float* Matrix::operator[](unsigned r)
{
    return &Values[r*NC];
}

Matrix::~Matrix()
{
    delete[] Values;
}
#endif

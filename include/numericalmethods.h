#ifndef NUMERICALMETHODS_H
#define NUMERICALMETHODS_H

//#include <QObject>
#include <vector>
#include <array>
#include <iostream>
#include "../include/numericalmatrix.h"
#include <math.h>

template<class T>
class NumericalMethods
{
public:
    NumericalMethods();

    static double determinant(int n, NumericalMatrix<T> mat);
    static NumericalMatrix<T> inverse(NumericalMatrix<T> m_mat);
    static NumericalMatrix<T> transpose(NumericalMatrix<T> m_mat);


private:
    NumericalMatrix<T> m_mat;
    NumericalMatrix<T> m_vec;


};

template class NumericalMethods<int>;
template class NumericalMethods<double>;
#endif // NUMERICALMETHODS_H

#ifndef NUMERICALMETHODS_H
#define NUMERICALMETHODS_H

//#include <QObject>
#include <vector>
#include <array>
#include <iostream>
#include "../include/numericalmatrix.h"
#include <math.h>
#include <functional>

template<class T>
class NumericalMethods
{
public:
    NumericalMethods();
    ~NumericalMethods();

    static double determinant(int n, NumericalMatrix<T> mat);
    static NumericalMatrix<T> inverse(NumericalMatrix<T> mat);
    static NumericalMatrix<T> transpose(NumericalMatrix<T> mat);

    template<typename G>
    static NumericalMatrix<T> newtonsMethod(G&& f){
        // F[0] = f1; F[1] = f2;
        f.initializeFunctions(f.X, f.JJ,f.F,f.iJJ);
        // number of iterations
        int N = 100000;
        // current iteration
        int i = 0;

        // main loop to solve
        for (i = 0; i < N; i++){

            // jacobian with initial guess
            f.initializeJacobian(f.JJ);
            // inverse of jacobian with initial guess
            f.computeInverseJacobian(f.JJ,f.iJJ);
            // take original f1, f2 and compute guess.
            f.computeFunction(f.F);
            // compute new guessing points
            if(f.computeNewGuess(f.X,i)){
                i=N;
                return f.X;
            }
        }
    }

    template<typename F>
    static double trap(F&& f, double a, double b){

        int N = 10000;
        double step = (b-a)/N;
        double s = 0;
        for (int i=0; i<=N; i++) {
          double xi = a + i*step;
          if (i == 0 || i == N) { s += f(xi); }

          else { s += 2* f(xi); }
        }
        s *= (b-a)/(2*N);
        return s;

    }
private:
    NumericalMatrix<T> m_mat;
    NumericalMatrix<T> m_vec;


};

template class NumericalMethods<int>;
template class NumericalMethods<double>;
#endif // NUMERICALMETHODS_H

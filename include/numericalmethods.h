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

    static double determinant(int n, NumericalMatrix<T> mat);
    static NumericalMatrix<T> inverse(NumericalMatrix<T> m_mat);
    static NumericalMatrix<T> transpose(NumericalMatrix<T> m_mat);

//    struct Equations{

//        // need to define the set of functions
//        NumericalMatrix<T> F;
//        // next need to define the jacobians of F
//        NumericalMatrix<T> JJ;
//        // internal TODO: needs to be internal of solver
//        NumericalMatrix<T> iJJ;
//        // The guess vector
//        NumericalMatrix<T> X;
//        // error
//        double eps = .000001;

//        void initializeJacobian(NumericalMatrix<T> & JJ){
//            // jacobian with initial guess
//            JJ.setValue(0,0,2*X.getElement(0));
//            JJ.setValue(0,1,2*X.getElement(1));
//            JJ.setValue(1,0,2*X.getElement(1));
//            JJ.setValue(1,1,2*X.getElement(0));
//        }

//        void computeFunction(NumericalMatrix<T> & F){
//            F.setValue(0,pow(this->X.getElement(0), 2) + pow(this->X.getElement(1), 2) - 50);
//            F.setValue(1,this->X.getElement(0) * this->X.getElement(1) - 25);
//        }

//        void computeInverseJacobian(NumericalMatrix<T> JJ, NumericalMatrix<T> & iJJ){
//            iJJ = NumericalMethods<T>::inverse(JJ);

//        }

//        bool computeNewGuess(NumericalMatrix<T> &X, int i){
//            X.setValue(0,X.getElement(0) - (this->iJJ[0][0] * this->F.getElement(0) + this->iJJ[0][1] * this->F.getElement(1)));
//            X.setValue(1,X.getElement(1) - (this->iJJ[1][0] * this->F.getElement(1) + this->iJJ[1][1] * this->F.getElement(1)));

//            if ((abs(pow(X.getElement(0), 2) + pow(X.getElement(1), 2) - 50) < eps) && abs(X.getElement(0) * X.getElement(1) - 25)<eps)
//            {
//                printf("Solution converged at N = %3.3d to X[0] = %3.3f and X[1] = %3.3f\r\n",i, X.getElement(0), X.getElement(1));
////                return X;
//                return 1;
//            }
//        }
//    };

//    struct Equations;
    template<typename G>
    static NumericalMatrix<T> newtonsMethod(G&& f/*,NumericalMethods::Equations eqs*/){
        // F[0] = f1; F[1] = f2;
        NumericalMatrix<T> F = f.F;
        NumericalMatrix<T> JJ = f.JJ;
        NumericalMatrix<T> iJJ = f.iJJ;
        NumericalMatrix<T> X = f.X;

        // initial guess
        f.X.createVector({2,1});
        double determinant=0;

        // error
        double eps = .000001;
        // number of iterations
        int N = 100000;
        // current iteration
        int i = 0;

        f.JJ.createMatrix(2,2);
        f.F.createVector({0,0});
        f.iJJ.createMatrix(2,2);
        // main loop to solve
        for (i = 0; i < N; i++){

            // the sys of non eqs.
            //f1 = (x1)^2 +(x2)^2 - 50 =0
            //f2 = (x1)*(x2) -25=0

            // jacobian with initial guess
            f.initializeJacobian(f.JJ);
//            JJ.setValue(0,0,2*X.getElement(0));
//            JJ.setValue(0,1,2*X.getElement(1));
//            JJ.setValue(1,0,2*X.getElement(1));
//            JJ.setValue(1,1,2*X.getElement(0));

            determinant = NumericalMethods<T>::determinant(f.JJ.dimension(),f.JJ);

            // inverse of jacobian with initial guess
//            iJJ = NumericalMethods<T>::inverse(JJ);
            f.computeInverseJacobian(f.JJ,f.iJJ);



            // take original f1, f2 and compute guess.
            f.computeFunction(F);

//            F.setValue(0,pow(X.getElement(0), 2) + pow(X.getElement(1), 2) - 50);
//            F.setValue(1,X.getElement(0) * X.getElement(1) - 25);


            // compute new guessing points
            f.computeNewGuess(f.X,i);
//            X.setValue(0,X.getElement(0) - (iJJ[0][0] * F.getElement(0) + iJJ[0][1] * F.getElement(1)));
//            X.setValue(1,X.getElement(1) - (iJJ[1][0] * F.getElement(1) + iJJ[1][1] * F.getElement(1)));

//            if ((abs(pow(X.getElement(0), 2) + pow(X.getElement(1), 2) - 50) < eps) && abs(X.getElement(0) * X.getElement(1) - 25)<eps)
//            {
//                printf("Solution converged at N = %3.3i to X[0] = %3.3f and X[1] = %3.3f\r\n",i, X.getElement(0), X.getElement(1));
//                return X;
//            }
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

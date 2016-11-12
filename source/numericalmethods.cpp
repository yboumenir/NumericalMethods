//#include "../include/numericalmethods.h"

//template <class T>
//NumericalMethods<T>::NumericalMethods()
//{

//}
//template <class T>
//NumericalMethods<T>::~NumericalMethods()
//{

//}


//template <class T>
//NumericalMatrix<T> NumericalMethods<T>::inverse(NumericalMatrix<T> m_mat){
//}

//template <class T>
//NumericalMatrix<T> NumericalMethods<T>::transpose(NumericalMatrix<T> mat)
//{
//}
///*
//template <class T>
//NumericalMatrix<T> NumericalMethods<T>::newtonsMethod(NumericalMatrix<T> mat, NumericalMatrix<T> jacobs)
//{
//    // F[0] = f1; F[1] = f2;
//    NumericalMatrix<T> F;
//    NumericalMatrix<T> JJ;
//    NumericalMatrix<T> iJJ;
//    NumericalMatrix<T> X;

//    // initial guess
//    X.createVector({2,1});
//    double determinant=0;

//    // error
//    double eps = .000001;
//    // number of iterations
//    int N = 100000;
//    // current iteration
//    int i = 0;

//    JJ.createMatrix(2,2);
//    F.createVector({0,0});
//    iJJ.createMatrix(2,2);
//    // main loop to solve
//    for (i = 0; i < N; i++){

//        // the sys of non eqs.
//        //f1 = (x1)^2 +(x2)^2 - 50 =0
//        //f2 = (x1)*(x2) -25=0

//        // jacobian with initial guess
//        JJ.setValue(0,0,2*X.getElement(0));
//        JJ.setValue(0,1,2*X.getElement(1));
//        JJ.setValue(1,0,2*X.getElement(1));
//        JJ.setValue(1,1,2*X.getElement(0));

//        determinant = NumericalMethods<T>::determinant(JJ.dimension(),JJ);

//        // inverse of jacobian with initial guess
//        iJJ = NumericalMethods<T>::inverse(JJ);

//        // take original f1, f2 and compute guess.
//        F.setValue(0,pow(X.getElement(0), 2) + pow(X.getElement(1), 2) - 50);
//        F.setValue(1,X.getElement(0) * X.getElement(1) - 25);

//        // compute new guessing points
//        X.setValue(0,X.getElement(0) - (iJJ[0][0] * F.getElement(0) + iJJ[0][1] * F.getElement(1)));
//        X.setValue(1,X.getElement(1) - (iJJ[1][0] * F.getElement(1) + iJJ[1][1] * F.getElement(1)));

//        if ((abs(pow(X.getElement(0), 2) + pow(X.getElement(1), 2) - 50) < eps) && abs(X.getElement(0) * X.getElement(1) - 25)<eps)
//        {
//            printf("Solution converged at N = %3.3i to X[0] = %3.3f and X[1] = %3.3f\r\n",i, X.getElement(0), X.getElement(1));
//            return X;
//        }
//    }

//}
//*/

////// simple traps integration ...
////template <class T>
////template <typename F>
////NumericalMethods<T> NumericalMethods<T>::trap(F &&f, double a, double b)
////{
////      int N = 10000;
////      double step = (b-a)/N;
////      double s = 0;
////      for (int i=0; i<=N; i++) {
////        double xi = a + i*step;
////        if (i == 0 || i == N) { s += f(xi); }

////        else { s += 2* f(xi); }
////      }
////      s *= (b-a)/(2*N);
////      return s;
////}

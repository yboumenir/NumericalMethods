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
    NumericalMethods(){

    }

    ~NumericalMethods(){

    }

    static double determinant(int n, NumericalMatrix<T> mat){
        double d = 0;
        int c, subi, i, j, subj;
        NumericalMatrix<T> submat;
        submat.createMatrix(n,n);

        if (n == 2)
        {
            return((mat[0][0] * mat[1][1]) - (mat[1][0] * mat[0][1]));
        }
        else
        {
            for (c = 0; c < n; c++)
            {
                subi = 0;
                for (i = 1; i < n; i++)
                {
                    subj = 0;
                    for (j = 0; j < n; j++)
                    {
                        if (j == c)
                        {
                            continue;
                        }
                        submat.setValue(subi,subj,mat[i][j]);
                        subj++;
                    }
                    subi++;
                }
                d = d + (pow(-1, c) * mat[0][c] * determinant(n - 1, submat));
            }
        }
        return d;
}




    static NumericalMatrix<T> inverse(NumericalMatrix<T> mat){
        int i, j;

        double determinant = NumericalMethods<T>::determinant(mat.dimension(),mat);

        if (determinant == 0)
        {
            std::cout << "Inverse does not exist (Determinant=0).\n";
        }


        NumericalMatrix<T> A = mat;
        int n = A.dimension();
        NumericalMatrix<T> AInverse;
        AInverse.createMatrix(n,n);

        // A = input matrix (n x n)
        // n = dimension of A
        // AInverse = inverted matrix (n x n)
        // This function inverts a matrix based on the Gauss Jordan method.
        // The function returns 1 on success, 0 on failure.
        int iPass, imx, icol, irow;
        float det, temp, pivot, factor;
        det = 1;

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                AInverse.setValue(i,j,0);
            }
            AInverse.setValue(i,i,1);
            //expect identity here
    //        AInverse.printMatrix();
    //        AInverse[n*i+i] = 1;
        }

        // The current pivot row is iPass.
        // For each pass, first find the maximum element in the pivot column.
        for (iPass = 0; iPass < n; iPass++)
        {
            imx = iPass;
            for (irow = iPass; irow < n; irow++)
            {
                if (fabs(A[irow][iPass]) > fabs(A[imx][iPass]))
                    imx = irow;
            }
            // Interchange the elements of row iPass and row imx in both A and AInverse.
            if (imx != iPass)
            {
                for (icol = 0; icol < n; icol++)
                {
                    temp = AInverse[iPass][icol];
                    AInverse.setValue(iPass,icol,AInverse[imx][icol]);
                    AInverse.setValue(imx,icol,temp);

                    if (icol >= iPass)
                    {
                        temp = A[iPass][icol];
                        A.setValue(iPass,icol,A[imx][icol]);
                        A.setValue(imx,icol,temp);
                    }
                }
            }

            // The current pivot is now A[iPass][iPass].
            // The determinant is the product of the pivot elements.
            pivot = A[iPass][iPass];
            det = det * pivot;
            if (det == 0)
            {
                NumericalMatrix<T> ANULL;
                return (ANULL);
            }

            for (icol = 0; icol < n; icol++)
            {
                // Normalize the pivot row by dividing by the pivot element.
                AInverse.setValue(iPass,icol,(AInverse[iPass][icol] / pivot));
                if (icol >= iPass) A.setValue(iPass,icol,(A[iPass][icol] / pivot));
            }

            for (irow = 0; irow < n; irow++)
            // Add a multiple of the pivot row to each row.  The multiple factor
            // is chosen so that the element of A on the pivot column is 0.
            {
                if (irow != iPass) factor = A[irow][iPass];
                for (icol = 0; icol < n; icol++)
                {
                    if (irow != iPass)
                    {
                        AInverse.setValue(irow,icol,AInverse[irow][icol]- factor * AInverse[iPass][icol]);
                        A.setValue(irow,icol, (A[irow][icol]- factor * A[iPass][icol]));
                    }
                }
            }
        }

        return AInverse;

    }

    static NumericalMatrix<T> transpose(NumericalMatrix<T> mat){

        int size = mat.dimension();
        NumericalMatrix<T> tmp=mat;

        for(int i=0;i<size;i++){
            for (int j=0;j<size;j++){
            tmp.setValue(i,j,mat[j][i]);
            }
        }

    //    std::cout << "orig" << std::endl;
    //    m_mat.printMatrix();
    //    std::cout << "trans" << std::endl;
    //    tmp.printMatrix();
        return mat;
    }

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

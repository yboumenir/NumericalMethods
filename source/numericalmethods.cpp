#include "../include/numericalmethods.h"
template <class T>
NumericalMethods<T>::NumericalMethods()
{

}
template <class T>
NumericalMethods<T>::~NumericalMethods()
{

}

template <class T>
double NumericalMethods<T>::determinant(int n, NumericalMatrix<T>mat)
{
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

template <class T>
NumericalMatrix<T> NumericalMethods<T>::inverse(NumericalMatrix<T> m_mat){
    int i, j;

    double determinant = NumericalMethods<T>::determinant(m_mat.dimension(),m_mat);

    if (determinant == 0)
    {
        std::cout << "Inverse does not exist (Determinant=0).\n";
    }


    NumericalMatrix<T> A = m_mat;
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

template <class T>
NumericalMatrix<T> NumericalMethods<T>::transpose(NumericalMatrix<T> mat)
{
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
/*
template <class T>
NumericalMatrix<T> NumericalMethods<T>::newtonsMethod(NumericalMatrix<T> mat, NumericalMatrix<T> jacobs)
{
    // F[0] = f1; F[1] = f2;
    NumericalMatrix<T> F;
    NumericalMatrix<T> JJ;
    NumericalMatrix<T> iJJ;
    NumericalMatrix<T> X;

    // initial guess
    X.createVector({2,1});
    double determinant=0;

    // error
    double eps = .000001;
    // number of iterations
    int N = 100000;
    // current iteration
    int i = 0;

    JJ.createMatrix(2,2);
    F.createVector({0,0});
    iJJ.createMatrix(2,2);
    // main loop to solve
    for (i = 0; i < N; i++){

        // the sys of non eqs.
        //f1 = (x1)^2 +(x2)^2 - 50 =0
        //f2 = (x1)*(x2) -25=0

        // jacobian with initial guess
        JJ.setValue(0,0,2*X.getElement(0));
        JJ.setValue(0,1,2*X.getElement(1));
        JJ.setValue(1,0,2*X.getElement(1));
        JJ.setValue(1,1,2*X.getElement(0));

        determinant = NumericalMethods<T>::determinant(JJ.dimension(),JJ);

        // inverse of jacobian with initial guess
        iJJ = NumericalMethods<T>::inverse(JJ);

        // take original f1, f2 and compute guess.
        F.setValue(0,pow(X.getElement(0), 2) + pow(X.getElement(1), 2) - 50);
        F.setValue(1,X.getElement(0) * X.getElement(1) - 25);

        // compute new guessing points
        X.setValue(0,X.getElement(0) - (iJJ[0][0] * F.getElement(0) + iJJ[0][1] * F.getElement(1)));
        X.setValue(1,X.getElement(1) - (iJJ[1][0] * F.getElement(1) + iJJ[1][1] * F.getElement(1)));

        if ((abs(pow(X.getElement(0), 2) + pow(X.getElement(1), 2) - 50) < eps) && abs(X.getElement(0) * X.getElement(1) - 25)<eps)
        {
            printf("Solution converged at N = %3.3i to X[0] = %3.3f and X[1] = %3.3f\r\n",i, X.getElement(0), X.getElement(1));
            return X;
        }
    }

}
*/

//// simple traps integration ...
//template <class T>
//template <typename F>
//NumericalMethods<T> NumericalMethods<T>::trap(F &&f, double a, double b)
//{
//      int N = 10000;
//      double step = (b-a)/N;
//      double s = 0;
//      for (int i=0; i<=N; i++) {
//        double xi = a + i*step;
//        if (i == 0 || i == N) { s += f(xi); }

//        else { s += 2* f(xi); }
//      }
//      s *= (b-a)/(2*N);
//      return s;
//}

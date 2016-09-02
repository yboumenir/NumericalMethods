#include "../include/numericalmethods.h"
template <class T>
NumericalMethods<T>::NumericalMethods()
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
//    /*float* A, int n,*/ float* AInverse;
    NumericalMatrix<T> AInverse;
    AInverse.createMatrix(n,n);

    // A = input matrix (n x n)
    // n = dimension of A
    // AInverse = inverted matrix (n x n)
    // This function inverts a matrix based on the Gauss Jordan method.
    // The function returns 1 on success, 0 on failure.
    int iPass, imx, icol, irow;
    float det, temp, pivot, factor;
//    float* ac = (float*)calloc(n*n, sizeof(float));
    det = 1;

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            AInverse.setValue(i,j,0);
//            AInverse[n*i+j] = 0;
//            ac[n*i+j] = A[n*i+j];
//            ac[n*i+j] = A[i][j];
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
//                temp = AInverse[n*iPass+icol];
                temp = AInverse[iPass][icol];
//                AInverse[n*iPass+icol] = AInverse[n*imx+icol];
                AInverse.setValue(iPass,icol,AInverse[imx][icol]);

//                AInverse[n*imx+icol] = temp;
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
//            free(ac);
            NumericalMatrix<T> ANULL;
            return (ANULL);
        }

        for (icol = 0; icol < n; icol++)
        {
            // Normalize the pivot row by dividing by the pivot element.
//            AInverse[n*iPass+icol] = AInverse[n*iPass+icol] / pivot;
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
//                    AInverse[n*irow+icol] -= factor * AInverse[n*iPass+icol];
                    AInverse.setValue(irow,icol,AInverse[irow][icol]- factor * AInverse[iPass][icol]);
                    A.setValue(irow,icol, (A[irow][icol]- factor * A[iPass][icol]));
                }
            }
        }
    }

//    AInverse.printMatrix();

//    free(ac);
    return AInverse;
//    return 1;
//    for (j = 0; j<3; j++)
//    {
//        cout << ((a[(i + 1) % 3][(j + 1) % 3] *a[(i + 2) % 3][(j + 2) % 3]) - (a[(i + 1) % 3][(j + 2) % 3] *	a[(i + 2) % 3][(j + 1) % 3])) / determinant << "\t";
//    }
//    cout << "\n";


}

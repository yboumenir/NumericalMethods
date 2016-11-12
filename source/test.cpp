#include "../include/numericalmatrix.h"
#include "../include/numericalmethods.h"

struct Functions{

    static double foo(double x){
        return pow(cos(x),2)+sin(x);
    }

    // need to define the set of functions
    NumericalMatrix<double> F;
    // next need to define the jacobians of F
    NumericalMatrix<double> JJ;
    // internal TODO: needs to be internal of solver
    NumericalMatrix<double> iJJ;
    // The guess vector
    NumericalMatrix<double> X;
    void initializeFunctions(NumericalMatrix<double> & X,
                             NumericalMatrix<double>& JJ,
                             NumericalMatrix<double> & F,
                             NumericalMatrix<double> & iJJ){
        X.createVector({2,1});
        JJ.createMatrix(2,2);
        F.createVector({0,0});
        iJJ.createMatrix(2,2);
    }

    // error
    double eps = .000001;

    void initializeJacobian(NumericalMatrix<double> & JJ){
        // jacobian with initial guess
        JJ.setValue(0,0,2*X.getElement(0));
        JJ.setValue(0,1,2*X.getElement(1));
        JJ.setValue(1,0,2*X.getElement(1));
        JJ.setValue(1,1,2*X.getElement(0));
    }
    // the sys of non eqs.
    //f1 = (x1)^2 +(x2)^2 - 50 =0
    //f2 = (x1)*(x2) -25=0
    void computeFunction(NumericalMatrix<double> & F){
        this->F.setValue(0,pow(this->X.getElement(0), 2) + pow(this->X.getElement(1), 2) - 50);
        this->F.setValue(1,this->X.getElement(0) * this->X.getElement(1) - 25);
    }

    void computeInverseJacobian(NumericalMatrix<double> JJ, NumericalMatrix<double> & iJJ){
        iJJ = NumericalMethods<double>::inverse(JJ);

    }

    bool computeNewGuess(NumericalMatrix<double> &X, int i){
        X.setValue(0,X.getElement(0) - (this->iJJ[0][0] * this->F.getElement(0) + this->iJJ[0][1] * this->F.getElement(1)));
        X.setValue(1,X.getElement(1) - (this->iJJ[1][0] * this->F.getElement(1) + this->iJJ[1][1] * this->F.getElement(1)));

        if ((abs(pow(X.getElement(0), 2) + pow(X.getElement(1), 2) - 50) < eps) && abs(X.getElement(0) * X.getElement(1) - 25)<eps)
        {
            printf("Solution converged at N = %3.3i to X[0] = %3.3f and X[1] = %3.3f\r\n",i, X.getElement(0), X.getElement(1));
            return 1;
        }
    }


} funcs;

#include <functional>
#include <memory>
int main()
{
    try {
        NumericalMatrix<double>  mat;
        int row = 3;
        int col = 3;

        mat.createMatrix(row,col);

        mat.pushColumns(0,{7,2,1});
        mat.pushColumns(1,{0,3,-1});
        mat.pushColumns(2,{-3,4,-2});

        std::cout << "A: " << std::endl;
        mat.printMatrix();

        NumericalMethods<double>::transpose(mat);

        std::cout << "det: " << NumericalMethods<double>::determinant(mat.dimension(),mat) << std::endl;

        std::cout << "inv: "<< std::endl;
        NumericalMethods<double>::inverse(mat).printMatrix();

        std::cout << "newton solver" << std::endl;
        NumericalMethods<double>::newtonsMethod(funcs);

        std::cout << "integrating cos(x) from 0 to 3.14/2  :" << std::endl;
        std::cout << NumericalMethods<double>::trap(Functions::foo,0,M_PI/2);
        std::cout << std::endl;

//        intStack.createVector({1,2});

//        // manipulate int stack
//        intStack.push(7);
//        intStack.printVector();
//        std::cout << intStack.top() <<std::endl;


    }
    catch (std::exception const& ex) {
        std::cerr << "Exception: " << ex.what() <<std::endl;
        return -1;
    }
}


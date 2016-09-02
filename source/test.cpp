#include <iostream>
#include "../include/numericalmatrix.h"
#include "../include/numericalmethods.h"
int main()
{
    try {
        NumericalMatrix<int>       intStack;  // stack of ints
        NumericalMatrix<double>  mat;
        int row = 3;
        int col = 3;

        mat.createMatrix(row,col);
//        for(int i =0 ; i < row; i++ ){
//            for(int j=0; j <col; j++){
//                mat.setValue(j,i,rand()/1000000);
//            }
//        }

        mat.pushColumns(0,{7,2,1});
        mat.pushColumns(1,{0,3,-1});
        mat.pushColumns(2,{-3,4,-2});

        std::cout << "A: " << std::endl;
        mat.printMatrix();

        std::cout << "det: " << NumericalMethods<double>::determinant(mat.dimension(),mat) << std::endl;

        std::cout << "inv: "<< std::endl;
        NumericalMethods<double>::inverse(mat).printMatrix();


        intStack.createVector({1,2});

        // manipulate int stack
        intStack.push(7);
        intStack.printVector();
        std::cout << intStack.top() <<std::endl;

        // manipulate string stack
//        stringStack.push("hello");
//        std::cout << stringStack.top() << std::endl;
//        stringStack.pop();
//        stringStack.pop();
    }
    catch (std::exception const& ex) {
        std::cerr << "Exception: " << ex.what() <<std::endl;
        return -1;
    }
}


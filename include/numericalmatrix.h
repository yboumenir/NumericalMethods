#ifndef NUMERICALMATRIX_H
#define NUMERICALMATRIX_H

#include <vector>
#include <iostream>
#include <forward_list>
#include <stdexcept>

template<class T>
class NumericalMatrix
{
private:
    std::vector<std::vector<T>> m_matrix;
    std::vector<T>         m_vector;

public:
    NumericalMatrix(){
        m_matrix;
        m_vector;
    }

    ~NumericalMatrix(){

    }

    void pushColumns(int row,std::initializer_list<T> vec){
        for(auto v:vec){
            m_matrix[row]=vec;
        }
    }

    void createMatrix(int rows,int cols){
        std::vector<std::vector<T>> matrix(rows, std::vector<T>(cols));
        m_matrix=matrix;
    }

    void createVector(std::initializer_list<T> vec){
        m_vector.insert(m_vector.end(),vec.begin(),vec.end());
    }

    void printMatrix(){
        int tmp_counter=0;
        int rows = m_matrix.size();
    //    std::cout << "rows:" << rows << std::endl;
    //    std::cout << "cols:" << m_matrix[0].size() << std::endl;
        for (int i=0;i<rows;i++){
            std::cout << "[" ;
            tmp_counter=0;
            for(auto s:m_matrix[i]){
                tmp_counter++;
                std::cout << s;
                if(tmp_counter!=m_matrix[i].size())
                    std::cout << ",";
            }
            std::cout << "]" << std::endl;
        }
    }

    void printVector(){
        int tmp_counter=0;
        int cols = m_vector.size();
        std::cout << "columns:" << cols << std::endl << "[";

        for(auto s:m_vector){
            tmp_counter++;
            std::cout << s;
            if(tmp_counter!=m_vector.size())
                std::cout << ",";
            else{
            std::cout << "]" << std::endl;
            }
        }
    }

    void push(T const& elem){
         m_vector.push_back(elem);
    }
        // push element
    void pop(){
        if (m_matrix.empty()) {
            throw std::out_of_range("Stack<>::pop(): empty stack");
        }
        // remove last element
        m_matrix.pop_back();
    }               // pop element
    T top() const{
        if (m_matrix.empty()) {
            throw std::out_of_range("Stack<>::top(): empty stack");
        }
        // return copy of last element
        return m_matrix[0].back();
    }
        // return top element
    bool empty() const{       // return true if empty.
        return m_matrix.empty();
    }
    std::vector<std::vector<T>> matrix() const
    {
        return m_matrix;
    }


    class Proxy {
    public:
        Proxy(std::vector<T> _array) : _array(_array) { }

        T operator[](int col) {
            return _array[col];
        }
    private:
        std::vector<T> _array;
    };

    Proxy operator[](int row) {
        return Proxy(m_matrix[row]);
    }


    void setValue(int row, int col, T value){
        m_matrix[row][col] = value;
    }

    void setValue(int row, T value){
        m_vector[row] = value;
    }

    T    getElement(int row){
        return m_vector[row];
    }

    double dimension(){
        return matrix().size();
    }

};



// template grbg
//template class NumericalMatrix<int>;
//template class NumericalMatrix<double>;
//template class NumericalMatrix<float>;
//template class NumericalMatrix<std::string>;


#endif // NUMERICALMATRIX_H

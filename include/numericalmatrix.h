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
  NumericalMatrix();
  ~NumericalMatrix();
  void pushColumns(int row,std::initializer_list<T>);
  void createMatrix(int rows,int cols);
  void createVector(std::initializer_list<T>);
  void printMatrix();
  void printVector();
  void push(T const&);  // push element
  void pop();               // pop element
  T top() const;            // return top element
  bool empty() const{       // return true if empty.
      return m_matrix.empty();
  }
  std::vector<std::vector<T> > matrix() const;

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


  void setValue(int row, int col, T value);
  void setValue(int row, T value);
  T    getElement(int row);
  double dimension();

};



// template grbg
template class NumericalMatrix<int>;
template class NumericalMatrix<double>;
template class NumericalMatrix<float>;
//template class NumericalMatrix<std::string>;


#endif // NUMERICALMATRIX_H

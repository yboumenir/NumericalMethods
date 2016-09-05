#include "../include/numericalmatrix.h"
template <class T>
std::vector<std::vector<T>> NumericalMatrix<T>::matrix() const
{
    return m_matrix;
}
template <class T>
void NumericalMatrix<T>::setValue(int row, int col, T value)
{
    m_matrix[row][col] = value;

}
template <class T>
void NumericalMatrix<T>::setValue(int row, T value)
{
    m_vector[row] = value;
}
template <class T>
T NumericalMatrix<T>::getElement(int row)
{
    return m_vector[row];
}

template <class T>
double NumericalMatrix<T>::dimension()
{
    return matrix().size();
}
template <class T>
NumericalMatrix<T>::NumericalMatrix()
{
    m_matrix;
    m_vector;
}
template <class T>
NumericalMatrix<T>::~NumericalMatrix()
{

}
template <class T>
void NumericalMatrix<T>::pushColumns(int row, std::initializer_list<T> vec)
{
    for(auto v:vec){
        m_matrix[row]=vec;
    }
//    m_matrix.insert(m_vector.end(),vec.begin(),vec.end());
}
template <class T>
void NumericalMatrix<T>::createMatrix(int rows, int cols)
{
    std::vector<std::vector<T>> matrix(rows, std::vector<T>(cols));
    m_matrix=matrix;
}

template <class T>
void NumericalMatrix<T>::createVector(std::initializer_list<T> vec)
{
    m_vector.insert(m_vector.end(),vec.begin(),vec.end());
}

template <class T>
void NumericalMatrix<T>::printMatrix()
{
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
template <class T>
void NumericalMatrix<T>::printVector()
{
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

template <class T>
void NumericalMatrix<T>::push (T const& elem)
{
    // append copy of passed element
    m_vector.push_back(elem);

}

template <class T>
void NumericalMatrix<T>::pop ()
{
    if (m_matrix.empty()) {
        throw std::out_of_range("Stack<>::pop(): empty stack");
    }
    // remove last element
    m_matrix.pop_back();
}

template <class T>
T NumericalMatrix<T>::top () const
{
    if (m_matrix.empty()) {
        throw std::out_of_range("Stack<>::top(): empty stack");
    }
    // return copy of last element
    return m_matrix[0].back();
}

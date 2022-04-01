#ifndef _MATRIX_H_
#define _MATRIX_H_
#include <cassert>
#include <cstddef>
#include <cstdio>


template <typename T> 
class CMatrix {
public:
    CMatrix();
    /**
     * @brief 析构函数
     *
     */
    ~CMatrix();

public:
    void Resize(unsigned int row, unsigned int col);
    /**
     * @brief 保证矩阵不超过规定行列数
     *
     * @param row
     * @param col
     */
    inline T &operator()(unsigned int row, unsigned int col) {
      
        assert(row < m_Row && col < m_Col);
        return m_A[row][col];
    }

public:
     CMatrix& operator+(CMatrix& v)
    {
        for (int i=0;i!=m_Row;i++)
        {
            for (int j=0;j!=m_Col;j++)
            {
                m_A[i][j]=v(i,j)+m_A[i][j];
            }
        }
        return (*this);
    } 
    CMatrix& operator=(CMatrix& v)
    {
        for (int i=0;i!=m_Row;i++)
        {
            for (int j=0;j!=m_Col;j++)
            {
                m_A[i][j]=v(i,j);
            }
        }
        return (*this);
    }
    //  CMatrix& operator=(CMatrix&& v)
    // {
    //     for (int i=0;i!=m_Row;i++)
    //     {
    //         for (int j=0;j!=m_Col;j++)
    //         {
    //             m_A[i][j]=v(i,j);
    //         }
    //     }
    //     return (*this);
    // }
    template <typename T1>
    friend CMatrix<T1>& operator*(T1 x,CMatrix<T1>& c_temp);

//private:
    T **m_A;
    unsigned int m_Row;
    unsigned int m_Col;
};
template <typename T>
CMatrix<T>& operator*(T x,CMatrix<T>& c_temp)
{
    for (int i=0;i!=c_temp.m_Row;i++)
    {
        for (int j=0;j!=c_temp.m_Col;j++)
        {
            c_temp(i,j)=x*c_temp(i,j);
        }
    }
    return c_temp;
}  

///////////////////////////////////////////////////////////
//////////////
template <typename T>
CMatrix<T>::CMatrix(void) : m_Row(0), m_Col(0), m_A(NULL) {}
/**
 * @brief 析构函数
 *
 * @tparam T
 */
template <typename T> CMatrix<T>::~CMatrix(void) {
    for (unsigned int i = 0; i != m_Row; ++i)
        delete[] m_A[i];

    delete[] m_A;

    m_A = NULL;
}
/**
 * @brief 给矩阵分配内存
 *
 * @tparam T
 * @param row
 * @param col
 */
template <typename T>
void CMatrix<T>::Resize(unsigned int row, unsigned int col) {
    if (row == 0 || col == 0)
        return;

    m_Row = row;
    m_Col = col;
    m_A = new T *[m_Row];
    for (unsigned int i = 0; i != m_Row; ++i)
        m_A[i] = new T[m_Col];
}
#endif
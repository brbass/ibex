#ifndef Matrix_Functions_hh
#define Matrix_Functions_hh

#include "Check.hh"

#include <vector>

/*
  Templated functions that act on matrices
  All functions use column-first indexing
  The functions return vectors directly to avoid issues with overwriting
  memory when given and return vectors are the same
*/
namespace Matrix_Functions
{
    template<class T> vector<T> transpose(int i_size, // rows
                                          int j_size, // cols
                                          vector<T> const &a)
    {
        Check(a.size() == i_size * j_size);
        
        vector<T> x(i_size * j_size);
        for (int i = 0; i < i_size; ++i)
        {
            for (int j = 0; j < j_size; ++j)
            {
                x[i + i_size * j] = a[j + j_size * i];
            }
        }
        
        return x;
    }
    
    template<class T> vector<T> square_matrix_vector_product(int size,
                                                             vector<T> const &a,
                                                             vector<T> const &b)
    {
        Check(a.size() = size * size);
        Check(b.size() == size);
        
        vector<T> x(size);
        for (int i = 0; i < size; ++i)
        {
            T sum = b[0] * a[size * i];
            for (int j = 1; j < size; ++j)
            {
                sum += a[j + size * i] * b[j];
            }
            x[i] = sum;
        }

        return x;
    }

    template<class T> vector<T> square_matrix_matrix_product(int size,
                                                             vector<T> const &a,
                                                             vector<T> const &b)
    {
        Check(a.size() == size * size);
        Check(b.size() == size * size);

        vector<double> x(size * size);
        for (int i = 0; i < size; ++i)
        {
            for (int j = 0; j < size; ++j)
            {
                T sum = a[size * i] * b[j];
                for (int k = 1; k < size; ++k)
                {
                    sum += a[k + size * i] * b[j + size * k];
                }
                x[j + size * i] = sum;
            }
        }
        
        return x;
    }

    template<class T> vector<T> matrix_matrix_product(int i_size, // a rows
                                                      int k_size, // a cols / b rows
                                                      int j_size, // b cols
                                                      vector<T> const &a,
                                                      vector<T> const &b);
    {
        Check(a.size() == i_size * k_size);
        Check(b.size() == k_size * j_size);
        
        vector<T> x(i_size * j_size);
        for (int i = 0; i < i_size; ++i)
        {
            for (int j = 0; j < j_size; ++j)
            {
                T sum = a[k_size * i] * b[j];
                for (int k = 1; k < k_size; ++k)
                {
                    sum += a[k + k_size * i] * b[j + j_size * k];
                }
                x[j + j_size * i] = sum;
            }
        }

        return x;
    }
} // end namespace Matrix_Functions

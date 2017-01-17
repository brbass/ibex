#ifndef Vector_Functions_hh
#define Vector_Functions_hh

#include "Check.hh"

#include <cmath>
#include <vector>

/*
  Templated functions that act on vectors
*/
namespace Vector_Functions
{
    // Vector/vector functions
    
    template<class T> std::vector<T> add(std::vector<T> const &x,
                                         std::vector<T> const &y)
    {
        int x_size = x.size();
        int y_size = y.size();
        
        Check(x_size == y_size);

        std::vector<T> result(x_size);
        
        for (int i = 0; i < x_size; ++i)
        {
            result[i] = x[i] + y[i];
        }
        
        return result;
    }

    template<class T> std::vector<T> subtract(std::vector<T> const &x,
                                              std::vector<T> const &y)
    {
        int x_size = x.size();
        int y_size = y.size();
        
        Check(x_size == y_size);

        std::vector<T> result(x_size);
        
        for (int i = 0; i < x_size; ++i)
        {
            result[i] = x[i] - y[i];
        }
        
        return result;
    }
    
    template<class T> T dot(std::vector<T> const &x,
                            std::vector<T> const &y)
    {
        int x_size = x.size();
        int y_size = y.size();
        
        Check(x_size == y_size);
        Check(x_size > 0);
        
        T result = x[0] * y[0];
        
        for (int i = 1; i < x_size; ++i)
        {
            result += x[i] * y[i];
        }
        
        return result;
    }

    template<class T> std::vector<T> tensor_dot(std::vector<T> const &x,
                                                std::vector<T> const &y)
    {
        int x_size = x.size();
        int y_size = y.size();

        Check(x_size % y_size == 0);
        Check(y_size > 0);

        int r_size = x_size / y_size;
        
        std::vector<T> result(r_size);
        
        for (int i = 0; i < r_size; ++i)
        {
            result[i] = x[0 + y_size * i] * y[0];
            
            for (int j = 1; j < y_size; ++j)
            {
                int k = j + y_size * i;
                
                result[i] += x[k] * y[j];
            }
        }
        
        return result;
    }
    
    template<class T> std::vector<T> tensor_product(std::vector<T> const &x,
                                                    std::vector<T> const &y)
    {
        int x_size = x.size();
        int y_size = y.size();

        Check(x_size == y_size);
        Check(x_size > 0);

        std::vector<T> result(x_size * x_size);

        for (int i = 0; i < x_size; ++i)
        {
            for (int j = 0; j < x_size; ++j)
            {
                int k = j + x_size * i;

                result[k] = x[i] * y[j];
            }
        }

        return result;
    }

    template<class T> std::vector<T> mean(std::vector<T> const &x,
                                          std::vector<T> const &y)
    {
        int x_size = x.size();
        int y_size = y.size();
        
        Check(x_size == y_size);
        Check(x_size > 0);
        
        std::vector<T> result(y_size);

        for (int i = 0; i < x_size; ++i)
        {
            result[i] = (x[i] + y[i]) / 2;
        }

        return result;
    }
    
    // Vector/scalar functions
    
    template<class T> std::vector<T> add(std::vector<T> const &x,
                                         T const y)
    {
        int x_size = x.size();
        
        std::vector<T> result(x_size);
        
        for (int i = 0; i < x_size; ++i)
        {
            result[i] = x[i] + y;
        }
        
        return result;
    }
    
    template<class T> std::vector<T> subtract(std::vector<T> const &x,
                                              T const y)
    {
        int x_size = x.size();
        
        std::vector<T> result(x_size);
        
        for (int i = 0; i < x_size; ++i)
        {
            result[i] = x[i] - y;
        }
        
        return result;
    }

    template<class T> std::vector<T> multiply(std::vector<T> const &x,
                                              T const t)
    {
        int x_size = x.size();
        
        std::vector<T> result(x_size);
        
        for (int i = 0; i < x_size; ++i)
        {
            result[i] = x[i] * t;
        }
        
        return result;
    }

    template<class T> std::vector<T> power(std::vector<T> const &x,
                                           T const t)
    {
        int x_size = x.size();
        
        std::vector<T> result(x_size);
        
        for (int i = 0; i < x_size; ++i)
        {
            result[i] = std::pow(x[i], t);
        }
        
        return result;
    }

    // Vector functions
    
    template<class T> std::vector<T> abs(std::vector<T> const &x)
    {
        int x_size = x.size();

        std::vector<T> result(x_size);

        for (int i = 0; i < x_size; ++i)
        {
            result[i] = std::abs(x[i]);
        }

        return result;
    }
    
    template<class T> T magnitude(std::vector<T> const &x)
    {
        int x_size = x.size();

        Check(x_size > 0);

        T result = x[0] * x[0];
        
        for (int i = 1; i < x_size; ++i)
        {
            result += x[i] * x[i];
        }
        
        return sqrt(result);
    }

    template<class T> T magnitude_squared(std::vector<T> const &x)
    {
        int x_size = x.size();

        Check(x_size > 0);

        T result = x[0] * x[0];
        
        for (int i = 1; i < x_size; ++i)
        {
            result += x[i] * x[i];
        }
        
        return result;
    }
    
    template<class T> std::vector<T> normalize(std::vector<T> const &x)
    {
        int x_size = x.size();

        Check(x_size > 0);

        T sum = x[0] * x[0];
        
        for (int i = 1; i < x_size; ++i)
        {
            sum += x[i] * x[i];
        }

        T norm = sqrt(sum);

        std::vector<T> result(x_size);
        
        for (int i = 0; i < x_size; ++i)
        {
            result[i] = x[i] / norm;
        }
        
        return result;
    }
    
    template<class T> T max(std::vector<T> const &x)
    {
        int x_size = x.size();

        Check(x_size > 0);
        
        T result = x[0];

        for (int i = 1; i < x_size; ++i)
        {
            if (x[i] > result)
            {
                result = x[i];
            }
        }

        return result;
    }

    template<class T> T min(std::vector<T> const &x)
    {
        int x_size = x.size();

        Check(x_size > 0);
        
        T result = x[0];

        for (int i = 1; i < x_size; ++i)
        {
            if (x[i] < result)
            {
                result = x[i];
            }
        }

        return result;
    }
} // namespace Vector_Functions

#endif

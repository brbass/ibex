#ifndef Unit_Test_hh
#define Unit_Test_hh

#include <cmath>
#include <limits>
#include <vector>

namespace Check_Equality
{
    template<class T> bool equal(T v1,
                                 T v2) 
    {
        return v1 == v2;
    }
    
    template<class T> bool approx(T v1,
                                  T v2,
                                  T tolerance) 
    {
        return std::abs(v1 - v2) <= tolerance;
    }
    
    template<class T> bool equal(std::vector<T> const &v1,
                                 std::vector<T> const &v2) 
    {
        return v1 == v2;
    }
    
    template<class T> bool approx(std::vector<T> const &v1,
                                  std::vector<T> const &v2,
                                  T tolerance) 
    {
        int size1 = v1.size();
        int size2 = v2.size();
        
        if (size1 != size2)
        {
            return false;
        }
        
        for (int i = 0; i < size1; ++i)
        {
            if (std::abs(v1[i] - v2[i]) > tolerance)
            {
                return false;
            }
        }
        
        return true;
    }

    template<class T> bool norm_approx(std::vector<T> const &v1,
                                       std::vector<T> const &v2,
                                       T tolerance,
                                       T &norm) 
    {
        int size1 = v1.size();
        int size2 = v2.size();
        
        if (size1 != size2)
        {
            return false;
        }

        double sum = 0;
        for (int i = 0; i < size1; ++i)
        {
            sum += std::abs(v1[i] - v2[i]);
        }
        sum /= size1;
        norm = sum;
        
        if (norm > tolerance)
        {
            return false;
        }
        else
        {
            return true;
        }
    }
}

#endif

#ifndef Vector_Functions_2D_hh
#define Vector_Functions_2D_hh

#include <cmath>
#include <vector>

namespace Vector_Functions_2D
{
    using namespace std;
    
    template<class T> vector<T> add(vector<T> const &x,
                                    vector<T> const &y)
    {
        return {x[0] + y[0],
                x[1] + y[1]};
    }
    
    template<class T> vector<T> subtract(vector<T> const &x,
                                         vector<T> const &y)
    {
        return {x[0] - y[0],
                x[1] - y[1]};
    }
    
    template<class T> vector<T> multiply(vector<T> const &x,
                                         T const t)
    {
        return {x[0] * t,
                x[1] * t};
    }
    
    template<class T> T dot(vector<T> const &x,
                            vector<T> const &y)
    {
        return x[0] * y[0] + x[1] * y[1];
    }
    
    template<class T> T cross(vector<T> const &x,
                              vector<T> const &y)
    {
        return x[0] * y[1] - x[1] * y[0];
    }
    
    template<class T> T magnitude(vector<T> const &x)
    {
        return sqrt(x[0] * x[0] + x[1] * x[1]);
    }
    
    template<class T> T magnitude_squared(vector<T> const &x)
    {
        return x[0] * x[0] + x[1] * x[1];
    }

    template<class T> vector<T> normalize(vector<T> const &x)
    {
        double k = sqrt(x[0] * x[0] + x[1] * x[1]);
        
        return {x[0] / k,
                x[1] / k};
    }
} // namespace Vector_Functions_2D

#endif

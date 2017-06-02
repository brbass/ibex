#ifndef Conversion_hh
#define Conversion_hh

#include <string>
#include <utility>
#include <vector>

#include "Check.hh"

template<class T1, class T2>
class Conversion
{
public:
    
    Conversion(std::vector<std::pair<T1, T2> > const &conversions):
        conversions_(conversions)
    {
    }

    T1 convert(T2 value) const
    {
        for (std::pair<T1, T2> const &conversion : conversions_)
        {
            if (conversion.second == value)
            {
                return conversion.first;
            }
        }
        
        AssertMsg(false, "enum value (" + value + ") not found");
        return conversions_[0].first;
    }
    T2 convert(T1 value)
    {
        for (std::pair<T1, T2> const &conversion : conversions_)
        {
            if (conversion.first == value)
            {
                return conversion.second;
            }
        }

        AssertMsg(false, "enum value not found");
        return conversions_[0].second;
    }
    T1 convert(T2 value,
               T1 default_value) const
    {
        for (std::pair<T1, T2> const &conversion : conversions_)
        {
            if (conversion.second == value)
            {
                return conversion.first;
            }
        }
        
        return default_value;
    }
    T2 convert(T1 value,
               T2 default_value)
    {
        for (std::pair<T1, T2> const &conversion : conversions_)
        {
            if (conversion.first == value)
            {
                return conversion.second;
            }
        }

        return default_value;
    }
    
private:

    std::vector<std::pair<T1, T2> > conversions_;
};

#endif

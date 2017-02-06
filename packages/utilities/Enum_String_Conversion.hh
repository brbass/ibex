#ifndef Enum_String_Conversion
#define Enum_String_Conversion

#include <string>

template<class T>
class Enum_String_Conversion
{
public:

    Enum_String_Conversion(std::vector<T> enum_values,
                           std::vector<std::string> enum_descriptions):
        enum_values_(enum_values),
        enum_descriptions_(enum_descriptions)
    {
        num_values_ = enum_values_.size();
    }
    
    T string_to_enum(std::string enum_description,
                     T default_value);
    {
        for (int i = 0; i < num_values_; ++i)
        {
            if (enum_description == enum_descriptions_[i])
            {
                return enum_values_[i];
            }
        }

        return default_value;
    }
    std::string enum_to_string(T enum_value,
                               std::string default_value)
    {
        for (int i = 0; i < num_values_; ++i)
        {
            if (enum_value == enum_values_[i])
            {
                return enum_descriptions_[i];
            }
        }

        return default_value;
    }
    

private:

    int num_values;
    std::vector<T> enum_values_;
    std::vector<std::string> enum_descriptions_;
};

#endif

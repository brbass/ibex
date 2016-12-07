#ifndef XML_Functions_hh
#define XML_Functions_hh

#include <string>
#include <vector>

#include "pugixml.hh"

#include "String_Functions.hh"

namespace XML_Functions
{
    // Declarations

    template<typename T> T attr_value(pugi::xml_attribute attr);
    template<typename T> T text_value(pugi::xml_text text);
    template<typename T> vector<T> text_vector(pugi::xml_text text);
    
    // Specializations
    template<> bool attr_value<bool>(pugi::xml_attribute attr)
    {
        return attr.as_bool();
    }
    template<> int attr_value<int>(pugi::xml_attribute attr)
    {
        return attr.as_int();
    }
    template<> unsigned attr_value<unsigned>(pugi::xml_attribute attr)
    {
        return attr.as_uint();
    }
    template<> double attr_value<double>(pugi::xml_attribute attr)
    {
        return attr.as_double();
    }
    template<> string attr_value<string>(pugi::xml_attribute attr)
    {
        return static_cast<string>(attr.as_string());
    }
    template<> bool text_value<bool>(pugi::xml_text text)
    {
        return text.as_bool();
    }
    template<> int text_value<int>(pugi::xml_text text)
    {
        return text.as_int();
    }
    template<> unsigned text_value<unsigned>(pugi::xml_text text)
    {
        return text.as_uint();
    }
    template<> double text_value<double>(pugi::xml_text text)
    {
        return text.as_double();
    }
    template<> string text_value<string>(pugi::xml_text text)
    {
        return static_cast<string>(text.as_string());
    }
    
    template<typename T> vector<T> text_vector(pugi::xml_text text)
    {
        vector<T> value;
        
        String_Functions::string_to_vector(text_value<string>(text),
                                           value);
        
        return value;
    }
}

#endif

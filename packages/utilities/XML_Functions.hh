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
    template<typename T> std::vector<T> text_vector(pugi::xml_text text);

    // Definitions
    
    template<typename T> std::vector<T> text_vector(pugi::xml_text text)
    {
        std::vector<T> value;
        
        String_Functions::string_to_vector(text_value<std::string>(text),
                                           value);
        
        return value;
    }
    
    // Specializations
    
    template<> inline bool attr_value<bool>(pugi::xml_attribute attr)
    {
        return attr.as_bool();
    }
    template<> inline int attr_value<int>(pugi::xml_attribute attr)
    {
        return attr.as_int();
    }
    template<> inline unsigned attr_value<unsigned>(pugi::xml_attribute attr)
    {
        return attr.as_uint();
    }
    template<> inline double attr_value<double>(pugi::xml_attribute attr)
    {
        return attr.as_double();
    }
    template<> inline std::string attr_value<std::string>(pugi::xml_attribute attr)
    {
        return static_cast<std::string>(attr.as_string());
    }
    template<> inline bool text_value<bool>(pugi::xml_text text)
    {
        return text.as_bool();
    }
    template<> inline int text_value<int>(pugi::xml_text text)
    {
        return text.as_int();
    }
    template<> inline unsigned text_value<unsigned>(pugi::xml_text text)
    {
        return text.as_uint();
    }
    template<> inline double text_value<double>(pugi::xml_text text)
    {
        return text.as_double();
    }
    template<> inline std::string text_value<std::string>(pugi::xml_text text)
    {
        return static_cast<std::string>(text.as_string());
    }
    
}

#endif

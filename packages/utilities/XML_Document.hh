#ifndef XML_Document_hh
#define XML_Document_hh

#include <string>

#include "pugixml.hh"

#include "XML_Node.hh";

class XML_Document : public XML_Node
{
public:

    // Create empty XML document
    XML_Document();

    // Load XML document from file
    XML_Document(std::string name); 
    
private:
    
    pugi::xml_document xml_document_;
}

#endif

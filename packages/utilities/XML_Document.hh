#ifndef XML_Document_hh
#define XML_Document_hh

#include <string>

#include "XML_Node.hh"

class XML_Document : public XML_Node
{
public:

    // Create empty XML document
    XML_Document();

    // Load XML document from file
    XML_Document(std::string name);

    // Save document
    void save(std::string name);
    
private:

    // Handles the creation of an XML_Node
    XML_Document(std::shared_ptr<pugi::xml_document> xml_doc,
                 std::string name);

    // Data
    std::shared_ptr<pugi::xml_document> xml_doc_;
    
    
};

#endif

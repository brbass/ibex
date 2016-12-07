#include "XML_Document.hh"

#include "Check.hh"
#include "XML_Node.hh"

using namespace std;

namespace
{
    pugi::xml_node get_empty_document()
    {
        pugi::xml_document doc;

        return doc;
    }
    pugi::xml_node get_document(string filename)
    {
        pugi::xml_document doc;
        pugi::xml_parse_result res = doc.load_file(filename.c_str());
        
        if (!res)
        {
            string error_message
                = "Could not open xml input file \""
                + filename
                + "\"\terror message: "
                + res.description();
            AssertMsg(false, error_message);
        }
        
        return doc;
    }
}

XML_Document::
XML_Document():
    XML_Node(get_empty_document())
{
}

XML_Document::
XML_Document(string filename):
    XML_Node(get_document(filename))
{
}


#include "XML_Document.hh"

#include "pugixml.hh"

#include "Check.hh"
#include "XML_Node.hh"

using namespace std;

namespace
{
    shared_ptr<pugi::xml_document> get_empty_document()
    {
        shared_ptr<pugi::xml_document> doc = make_shared<pugi::xml_document>();

        return doc;
    }
    shared_ptr<pugi::xml_document> get_document(string filename)
    {
        shared_ptr<pugi::xml_document> doc = make_shared<pugi::xml_document>();
        pugi::xml_parse_result res = doc->load_file(filename.c_str());
        
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
    XML_Document(get_empty_document(),
                 "new_document")
{
}

XML_Document::
XML_Document(string filename):
    XML_Document(get_document(filename),
                 filename)
{
}

XML_Document::
XML_Document(shared_ptr<pugi::xml_document> xml_doc,
             string name):
    XML_Node(xml_doc,
             name),
    xml_doc_(xml_doc)
{
    
}

void XML_Document::
save(string name)
{
    xml_doc_->save_file(name.c_str());
}

string XML_Document::
path() const
{
    return xml_doc_->path();
}

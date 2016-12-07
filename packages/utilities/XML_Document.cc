#include "XML_Document.hh"

#include "Check.hh"
#include "XML_Node.hh"

using namespace std;

XML_Document::
XML_Document():
    XML_Node(xml_document_)
{
}

XML_Document::
XML_Document(string filename):
    XML_Document()
{
    if (!xml_document_.load_file(filename.c_str()))
    {
        AssertMsg(false, "Could not open xml input file \"" + filename + "\"");
    }
}


#include <iostream>
#include <stdexcept>
#include <string>

namespace ch_ns
{
    using namespace std;
    
    void check(std::string condition,
               std::string file,
               int line)
    {
        cerr << "condition \"";
        cerr << condition;
        cerr << "\" failed on line";
        cerr << line;
        cerr << " of \"";
        cerr << file;
        cerr << "\"";
        cerr << endl;
    
        throw runtime_error("check");
    }

    void check(std::string condition,
               std::string message,
               std::string file,
               int line)
    {
        cerr << "condition \"";
        cerr << condition;
        cerr << "\" failed with error \"";
        cerr << message;
        cerr << "\" on line ";
        cerr << line;
        cerr << " of \"";
        cerr << file;
        cerr << "\"";
        cerr << endl;
    
        throw runtime_error("check");
    }
} // namespace ch_ns

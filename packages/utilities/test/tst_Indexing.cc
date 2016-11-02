#include <iostream>

#include "Indexing.hh"

using namespace std;

int test_indexing()
{
    int checksum = 0;
    
    size_t const dimension = 6;
    vector<size_t> const size = {4, 302, 50, 7, 24, 2};
    vector<size_t> const initial_subscript = {2, 107, 32, 5, 14, 1};
    
    Indexing<size_t> indexing(dimension,
                              size);
    
    size_t index = indexing.subscript_to_index(initial_subscript);
    
    vector<size_t> subscript = indexing.index_to_subscript(index);
    
    for (size_t i = 0; i < dimension; ++i)
    {
        if (initial_subscript[i] != subscript[i])
        {
            cout << "subscript wrong in dimension ";
            cout << i;
            cout << endl;
            cout << "\tinitial and final values: ";
            cout << initial_subscript[i];
            cout << "\t";
            cout << subscript[i];
            cout << endl;
            checksum += 1;
        }
    }
    
    return checksum;
}

int main()
{
    int checksum = 0;

    checksum += test_indexing();
    
    return checksum;
}

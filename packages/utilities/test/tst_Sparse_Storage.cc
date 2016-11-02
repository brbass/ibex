#include <iostream>
#include <limits>
#include <stdexcept> 

#include "Check_Equality.hh"
#include "Sparse_Storage.hh"
#include "Symmetric_Sparse_Storage.hh"

using namespace std;
namespace ce = Check_Equality;

int test_sparse_storage()
{
    int checksum = 0;
    double tolerance = 10 * numeric_limits<double>::epsilon();
    
    int number_of_elements = 1000;

    Sparse_Storage<int, double> storage(number_of_elements);

    int index1 = 100;
    int index2 = 20;
    double value = 2.7;
    
    storage.add_element(index1, index2, value);

    // Check that element exists
    if (!storage.element_exists(index1, index2))
    {
        checksum += 1;
        cout << "sparse: element does not exist" << endl;
    }
    
    // Check that spurious elements don't exist
    if (storage.element_exists(index1, index2 - 1) || storage.element_exists(index1 + 1, index2))
    {
        checksum += 1;
        cout << "sparse: spurious element exists" << endl;
    }

    // Check value of element
    if (!ce::approx(storage.get_element(index1, index2), value, tolerance))
    {
        checksum += 1;
        cout << "sparse: value of element incorrect" << endl;
    }
    
    // Check for correct error
    try
    {
        storage.get_element(index1 + 1, index2);
        checksum += 1;
        cout << "sparse: error not thrown 1" << endl;
    }
    catch (const out_of_range &oor)
    {
    }

    try
    {
        storage.get_element(index1, index2 + 1);
        checksum += 1;
        cout << "sparse: error not thrown 2" << endl;
    }
    catch (const out_of_range &oor)
    {
    }
    
    try
    {
        storage.add_element(index1 + number_of_elements, index2, value);
        checksum += 1;
        cout << "sparse: error not thrown 3" << endl;
    }
    catch (const out_of_range &oor)
    {
    }
    
    return checksum;
}

int test_symmetric_sparse_storage()
{
    int checksum = 0;
    double tolerance = 10 * numeric_limits<double>::epsilon();
    
    int number_of_elements = 1000;

    Symmetric_Sparse_Storage<int, double> storage(number_of_elements);

    int index1 = 100;
    int index2 = 20;
    double value = 2.7;
    
    storage.add_element(index1, index2, value);

    // Check that element exists
    if (!storage.element_exists(index1, index2))
    {
        checksum += 1;
        cout << "symmetric: element does not exist" << endl;
    }
    if (!storage.element_exists(index2, index1))
    {
        checksum += 1;
        cout << "symmetric: transpose element does not exist" << endl;
    }

    // Check that spurious elements don't exist
    if (storage.element_exists(index1, index2 - 1) || storage.element_exists(index1 + 1, index2))
    {
        checksum += 1;
        cout << "symmetric: spurious element exists" << endl;
    }

    // Check value of element
    if (!ce::approx(storage.get_element(index1, index2), value, tolerance))
    {
        checksum += 1;
        cout << "symmetric: value of element incorrect" << endl;
    }
    if (!ce::approx(storage.get_element(index2, index1), value, tolerance))
    {
        checksum += 1;
        cout << "symmetric: value of transpose element incorrect" << endl;
    }

    // Check for correct error
    try
    {
        storage.get_element(index1 + 1, index2);
        checksum += 1;
        cout << "symmetric: error not thrown 1" << endl;
    }
    catch (const out_of_range &oor)
    {
    }

    try
    {
        storage.get_element(index1, index2 + 1);
        checksum += 1;
        cout << "symmetric: error not thrown 2" << endl;
    }
    catch (const out_of_range &oor)
    {
    }
    
    try
    {
        storage.add_element(index1 + number_of_elements, index2 + number_of_elements, value);
        checksum += 1;
        cout << "symmetric: error not thrown 3" << endl;
    }
    catch (const out_of_range &oor)
    {
    }

    return checksum;
}

int main()
{
    int checksum = 0;
    
    checksum += test_sparse_storage();
    checksum += test_symmetric_sparse_storage();
    
    return checksum;
}

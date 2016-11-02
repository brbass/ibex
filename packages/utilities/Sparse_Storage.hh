#ifndef Sparse_Storage_hh
#define Sparse_Storage_hh

#include <iomanip>
#include <iostream>
#include <unordered_map>
#include <vector>

// Stores columns in a vector of unordered_map
// Does range checking on row, but not column
// Throws out_of_range exception if row does not exist
template<class Ordinal, class Scalar>
class Sparse_Storage
{
public:
    
    // Constructors
    Sparse_Storage(Ordinal number_of_elements):
        number_of_elements_(number_of_elements)
    {
        data_.resize(number_of_elements);
    }

    virtual int number_of_elements() const
    {
        return number_of_elements_;
    }
    
    // Creates element if one does not exist
    // Overwrites element if one does exist
    virtual void add_element(Ordinal index1,
                             Ordinal index2,
                             Scalar value)
    {
        data_.at(index1)[index2] = value;
    }
    
    // Elements are unique, so this is either one or zero
    virtual bool element_exists(Ordinal index1,
                                Ordinal index2) const
    {
        if (data_.at(index1).count(index2) > 0)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    
    // Throws out_of_range exception if element does not exist
    virtual Scalar get_element(Ordinal index1,
                               Ordinal index2) const
    {
        return data_.at(index1).at(index2);
    }

    // Print matrix
    friend std::ostream& operator<< (std::ostream &os,
                                     const Sparse_Storage<Ordinal, Scalar> &ss)
    {
        using std::setw;
        using std::endl;
        
        Ordinal w = 16;
        os << setw(w) << "row";
        os << setw(w) << "column";
        os << setw(w) << "value";
        os << endl;
        for (Ordinal i = 0; i < ss.number_of_elements_; ++i)
        {
            for (std::pair<Ordinal, Scalar> data : ss.data_[i])
            {
                os << setw(w) << i;
                os << setw(w) << data.first;
                os << setw(w) << data.second;
                os << endl;
            }
        }
        
        return os;
    }

    
protected:
    
    // Number of rows
    Ordinal number_of_elements_;
    
    // Data
    std::vector<std::unordered_map<Ordinal, Scalar> > data_;
};

#endif

#ifndef Symmetric_Sparse_Storage_hh
#define Symmetric_Sparse_Storage_hh

#include <unordered_map>
#include <vector>

#include "Sparse_Storage.hh"

template<class Ordinal, class Scalar>
class Symmetric_Sparse_Storage : public Sparse_Storage<Ordinal, Scalar>
{
public:
    
    // Constructor
    Symmetric_Sparse_Storage(Ordinal number_of_elements):
        Sparse_Storage<Ordinal, Scalar>(number_of_elements)
    {
    }
    
    // Creates element if one does not exist
    // Overwrites element if one does exist
    virtual void add_element(Ordinal index1,
                             Ordinal index2,
                             Scalar value) override
    {
        int row, column;
        get_indices(index1,
                    index2,
                    row,
                    column);
        
        this->data_.at(row)[column] = value;
    }
    
    // Elements are unique, so this is either one or zero
    virtual bool element_exists(Ordinal index1,
                                Ordinal index2) const override
    {
        int row, column;
        get_indices(index1,
                    index2,
                    row,
                    column);

        if (this->data_.at(row).count(column) > 0)
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
                               Ordinal index2) const override
    {
        int row, column;
        get_indices(index1,
                    index2,
                    row,
                    column);

        return this->data_.at(row).at(column);
    }
    
private:
    
    virtual void get_indices(Ordinal index1,
                             Ordinal index2,
                             Ordinal &row,
                             Ordinal &column) const
    {
        if (index1 < index2)
        {
            row = index1;
            column = index2;
        }
        else
        {
            row = index2;
            column = index1;
        }
    }
};

#endif

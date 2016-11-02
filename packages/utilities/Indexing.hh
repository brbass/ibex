#ifndef Indexing_hh
#define Indexing_hh

#include <cmath>
#include <vector>

template<class Ordinal>
class Indexing
{
public:
    
    Indexing(Ordinal dimension,
             std::vector<Ordinal> const &size):
        dimension_(dimension),
        size_(size)
    {
    }

    Ordinal subscript_to_index(std::vector<Ordinal> subscript)
    {
        Ordinal index = subscript[0];

        for (Ordinal d = 1; d < dimension_; ++d)
        {
            index = subscript[d] + size_[d] * index;
        }

        return index;
    }
    
    std::vector<Ordinal> index_to_subscript(Ordinal index)
    {
        std::vector<Ordinal> subscript(dimension_);

        Ordinal current_index = index;

        for (Ordinal i = dimension_ - 1; i > 0; --i)
        {
            subscript[i] = current_index % size_[i];
            
            current_index -= subscript[i];
            
            current_index /= size_[i];
        }
        
        subscript[0] = current_index;
        
        return subscript;
    }
    
private:
    
    Ordinal dimension_;
    std::vector<Ordinal> size_;
};

#endif

/* Old method
std::vector<Ordinal> index_to_subscript(Ordinal index)
{
    std::vector<Ordinal> subscript(dimension_);
        
    Ordinal product = 1;
    for (Ordinal i = 1; i < dimension_; ++i)
    {
        product *= size_[i];
    };
    
    Ordinal sum = index;
    for (Ordinal i = 0; i < dimension_ - 1; ++i)
    {
        subscript[i] = std::floor(static_cast<double>(sum) / product);
            
        sum -= product * subscript[i];
        product /= size_[i + 1];
    }
    
    subscript[dimension_ - 1] = sum;
        
    return subscript;
}
*/

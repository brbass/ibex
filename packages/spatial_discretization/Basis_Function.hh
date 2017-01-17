#ifndef Basis_Function_hh 
#define Basis_Function_hh

#include <memory>
#include <vector>

class Cartesian_Plane;
class Meshless_Function;

class Basis_Function
{
public:
    
    // Constructor
    Basis_Function(int index,
                   int dimension,
                   std::shared_ptr<Meshless_Function> meshless_function,
                   std::vector<std::shared_ptr<Cartesian_Plane> > boundary_surfaces);

    // Data access
    virtual int index() const
    {
        return index_;
    }
    virtual int dimension() const
    {
        return dimension_;
    }
    virtual int number_of_boundary_surfaces() const
    {
        return number_of_boundary_surfaces_;
    }
    virtual std::shared_ptr<Meshless_Function> function() const
    {
        return meshless_function_;
    }
    virtual std::shared_ptr<Cartesian_Plane> boundary_surface(int i) const
    {
        return boundary_surfaces_[i];
    }
    
private:
    
    // Data
    int index_;
    int dimension_;
    int number_of_boundary_surfaces_;
    std::shared_ptr<Meshless_Function> meshless_function_;
    std::vector<std::shared_ptr<Cartesian_Plane> > boundary_surfaces_;
};

#endif

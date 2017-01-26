#ifndef Boundary_Source_hh
#define Boundary_Source_hh

#include <memory>
#include <vector>

class Angular_Discretization;
class Energy_Discretization;
class XML_Node;

class Boundary_Source
{
public:

    struct Dependencies
    {
        enum class Angular
        {
            ISOTROPIC,
            MOMENTS,
            ORDINATES
        };

        Angular angular = Angular::ISOTROPIC;
    };
    
    Boundary_Source(int index,
                    Dependencies dependencies,
                    std::shared_ptr<Angular_Discretization> angular_discretization,
                    std::shared_ptr<Energy_Discretization> energy_discretization,
                    std::vector<double> const &boundary_source,
                    std::vector<double> const &alpha);
    
    int index() const
    {
        return index_;
    }
    int size() const
    {
        return size_;
    }
    Dependencies dependencies() const
    {
        return dependencies_;
    }
    std::vector<double> const &data() const
    {
        return boundary_source_;
    }
    std::vector<double> const &alpha() const
    {
        return alpha_;
    }
    std::shared_ptr<Angular_Discretization> angular_discretization() const
    {
        return angular_discretization_;
    }
    std::shared_ptr<Energy_Discretization> energy_discretization() const
    {
        return energy_discretization_;
    }
    bool has_reflection() const;
    void check_class_invariants() const;
    void output(XML_Node output_node) const;

private:

    int index_;
    int size_;
    Dependencies dependencies_;
    std::shared_ptr<Angular_Discretization> angular_discretization_;
    std::shared_ptr<Energy_Discretization> energy_discretization_;
    std::vector<double> boundary_source_;
    std::vector<double> alpha_;
};

#endif

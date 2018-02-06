#ifndef Heat_Transfer_Integration_hh
#define Heat_Transfer_Integration_hh

#include <memory>
#include <vector>

class Heat_Transfer_Data;
class Integration_Mesh;
class Integration_Mesh_Options;
class Weak_Spatial_Discretization;

struct Heat_Transfer_Integration_Options
{
    enum class Geometry
    {
        CYLINDRICAL_1D,
        CARTESIAN
    };

    Geometry geometry;
};

class Heat_Transfer_Integration
{
public:
    
    Heat_Transfer_Integration(std::shared_ptr<Heat_Transfer_Integration_Options> options,
                              std::shared_ptr<Integration_Mesh_Options> integration_options,
                              std::shared_ptr<Heat_Transfer_Data> data,
                              std::shared_ptr<Weak_Spatial_Discretization> spatial);

    // Data access
    std::vector<std::vector<double> > const &matrix() const
    {
        return matrix_;
    }
    std::vector<double> const &rhs() const
    {
        return rhs_;
    }
    
private:
    
    // Integration methods
    void initialize_integrals();
    void perform_integration();
    
    // Input data
    std::shared_ptr<Heat_Transfer_Integration_Options> options_;
    std::shared_ptr<Heat_Transfer_Data> data_;
    std::shared_ptr<Weak_Spatial_Discretization> spatial_;
    
    // Utility data
    std::shared_ptr<Integration_Mesh> mesh_;
    
    // Output data
    std::vector<std::vector<double> > matrix_;
    std::vector<double> rhs_;
};

#endif

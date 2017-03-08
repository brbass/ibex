#ifndef Weak_Spatial_Discretization_Factory_hh
#define Weak_Spatial_Discretization_Factory_hh

class Weak_Spatial_Discretization_Factory
{
public:

    Weak_Spatial_Discretization_Factory(std::shared_ptr<Solid_Geometry> solid_geometry,
                                        std::vector<std::shared_ptr<Cartesian_Plane> > const &boundary_surfaces);
    
    
    
private:

    vector<vector<double> > limits_;
    std::shared_ptr<Solid_Geometry> solid_geometry_;
    std::vector<std::shared_ptr<Cartesian_Plane> > boundary_surfaces_;
};

#endif

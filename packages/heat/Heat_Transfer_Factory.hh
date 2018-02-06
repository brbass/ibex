#ifndef Heat_Transfer_Factory_hh
#define Heat_Transfer_Factory_hh

#include <memory>
#include <string>

class Cartesian_Plane;
class Heat_Transfer_Data;
class Heat_Transfer_Solve;
class Solid_Geometry;
class Weak_Spatial_Discretization;

class Heat_Transfer_Factory
{
public:
    
    Heat_Transfer_Factory();

    void get_solid_1d(double length,
                      std::shared_ptr<Solid_Geometry> &solid,
                      std::vector<std::shared_ptr<Cartesian_Plane> > &boundary_surfaces) const;
    std::shared_ptr<Weak_Spatial_Discretization> get_spatial_discretization_1d(int number_of_points,
                                                                               double radius_num_intervals,
                                                                               double length,
                                                                               bool basis_mls,
                                                                               bool weight_mls,
                                                                               std::string basis_type,
                                                                               std::string weight_type) const;
    std::shared_ptr<Heat_Transfer_Solve> get_cylindrical_1d(int number_of_points,
                                                            double radius_num_intervals,
                                                            double length,
                                                            bool basis_mls,
                                                            bool weight_mls,
                                                            std::string basis_type,
                                                            std::string weight_type,
                                                            std::shared_ptr<Heat_Transfer_Data> data) const;
};

#endif

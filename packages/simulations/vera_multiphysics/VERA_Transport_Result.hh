#ifndef VERA_Transport_Result_hh
#define VERA_Transport_Result_hh

#include <memory>
#include <vector>

#include "Solver.hh"

class Angular_Discretization;
class Energy_Discretization;
class Solid_Geometry;
class Weak_Spatial_Discretization;
class XML_Node;

class VERA_Transport_Result
{
public:

    VERA_Transport_Result(int heat_dimension,
                          double fuel_radius,
                          double pincell_power,
                          std::shared_ptr<Solid_Geometry> solid,
                          std::shared_ptr<Angular_Discretization> angular,
                          std::shared_ptr<Energy_Discretization> energy,
                          std::shared_ptr<Weak_Spatial_Discretization> spatial,
                          std::shared_ptr<Solver> solver,
                          std::shared_ptr<Solver::Result> result);

    double get_fission_energy(std::vector<double> const &position);
    double get_radial_fission_energy(double radius);

    std::shared_ptr<Solver::Result> result() const
    {
        return result_;
    }
    
    void output_data(XML_Node output_node);
    
private:

    void normalize();

    // Input data
    int heat_dimension_;
    double pincell_power_;
    double fuel_radius_;
    std::shared_ptr<Solid_Geometry> solid_;
    std::shared_ptr<Angular_Discretization> angular_;
    std::shared_ptr<Energy_Discretization> energy_;
    std::shared_ptr<Weak_Spatial_Discretization> spatial_;
    std::shared_ptr<Solver> solver_;
    std::shared_ptr<Solver::Result> result_;
    int number_of_ordinates_;
    std::vector<double> ordinates_;
    std::vector<double> weights_;

    // Constants
    double const mev_to_joule_ = 1.6021766e-13;
    std::vector<double> const nu_ = {2.65063, 2.43223};
    std::vector<double> const kappa_ = {196.155, 193.083};
};

#endif

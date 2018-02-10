#ifndef VERA_Solid_Geometry_hh
#define VERA_Solid_Geometry_hh

#include <functional>
#include <memory>
#include <vector>

#include "Check.hh"
#include "Solid_Geometry.hh"
#include "XML_Node.hh"

class Angular_Discretization;
class Boundary_Source;
class Energy_Discretization;
class Material;
class Material_Factory;

class VERA_Solid_Geometry : public Solid_Geometry
{
public:

    VERA_Solid_Geometry(bool include_ifba,
                        std::function<double(std::vector<double> const &)> temperature,
                        std::shared_ptr<Angular_Discretization> angular,
                        std::shared_ptr<Energy_Discretization> energy,
                        std::vector<std::shared_ptr<Material> > materials,
                        std::shared_ptr<Boundary_Source> boundary_source);

    virtual int dimension() const override;
    virtual int find_region(std::vector<double> const &position) const override
    {
        AssertMsg(false, "not implemented");
        return Geometry_Errors::NO_REGION;
    }
    virtual int find_surface(std::vector<double> const &position) const override
    {
        AssertMsg(false, "not implemented");
        return Geometry_Errors::NO_SURFACE;
    }
    virtual int next_intersection(std::vector<double> const &initial_position,
                                  std::vector<double> const &initial_direction,
                                  int &final_region,
                                  double &distance,
                                  std::vector<double> &final_position) const override
    {
        AssertMsg(false, "not implemented");
        return Geometry_Errors::NO_SURFACE;
    }
    virtual int next_boundary(std::vector<double> const &initial_position,
                              std::vector<double> const &initial_direction,
                              int &boundary_region,
                              double &distance,
                              std::vector<double> &final_position) const override
    {
        AssertMsg(false, "not implemented");
        return Geometry_Errors::NO_SURFACE;
    }        
    virtual void optical_distance(std::vector<double> const &initial_position,
                                  std::vector<double> const &final_position,
                                  std::vector<double> &optical_distance) const override
    {
        AssertMsg(false, "not implemented");
    }
    virtual std::shared_ptr<Material> material(std::vector<double> const &position) const override;
    virtual std::shared_ptr<Boundary_Source> boundary_source(std::vector<double> const &position) const override;
    virtual void check_class_invariants() const override;
    virtual void output(XML_Node output_node) const override
    {
        AssertMsg(false, "not implemented");
    }
    
private:

    enum Material_Type
    {
        NONE = -1,
        FUEL = 0,
        IFBA = 1,
        GAP = 2,
        CLAD = 3,
        MOD = 4
    };
    enum Problem_Type
    {
        V1B = 0,
        V1E = 1
    };
    enum Temperature_Type
    {
        K600 = 0,
        K1000 = 1
    };
    
    double radial_distance2(std::vector<double> const &position) const;
    std::shared_ptr<Material> get_material_by_index(Material_Type mat_type,
                                                    Problem_Type prob_type,
                                                    Temperature_Type temp_type) const;
    std::vector<double> weighted_cross_section(double temperature,
                                               std::vector<double> const &xs600,
                                               std::vector<double> const &xs1000) const;
    std::shared_ptr<Material> weighted_material(double temperature,
                                                std::shared_ptr<Material> mat600,
                                                std::shared_ptr<Material> mat1000) const;
    
    // Data
    bool include_ifba_;
    std::function<double(std::vector<double> const &)> temperature_;
    std::shared_ptr<Angular_Discretization> angular_;
    std::shared_ptr<Energy_Discretization> energy_;
    std::shared_ptr<Material_Factory> material_factory_;
    std::vector<double> material_radii2_;

    // Cross sections
    int number_of_material_types_;
    int number_of_problem_types_;
    int number_of_temperature_types_;
    std::vector<std::shared_ptr<Material> > materials_;
    std::shared_ptr<Boundary_Source> boundary_source_;
};

#endif

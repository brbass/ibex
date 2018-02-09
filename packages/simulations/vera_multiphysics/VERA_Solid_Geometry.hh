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
                        std::shared_ptr<Energy_Discretization> energy);

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

    void initialize_materials();
    void get_cross_sections(double temperature,
                            std::vector<double> const &position,
                            std::vector<double> &sigma_t,
                            std::vector<double> &sigma_s,
                            std::vector<double> &sigma_f) const;

    // Data
    bool include_ifba_;
    std::function<double(std::vector<double> const &)> temperature_;
    std::shared_ptr<Angular_Discretization> angular_;
    std::shared_ptr<Energy_Discretization> energy_;
    std::shared_ptr<Material_Factory> material_factory_;

    // Cross sections
    enum Material_Type
    {
        FUEL = 0,
        IFBA = 1,
        GAP = 2,
        CLAD = 3,
        MOD = 4
    };
    enum XS_Type
    {
        NORM_600 = 0,
        NORM_1000 = 1,
        IFBA_600 = 2,
        IFBA_1000 = 3
    };
    int number_of_material_types_;
    int number_of_xs_types_;
    std::vector<std::vector<std::shared_ptr<Material> > > materials_;
};

#endif

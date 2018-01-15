#ifndef Manufactured_Solid_Geometry_hh
#define Manufactured_Solid_Geometry_hh

#include "Check.hh"
#include "Solid_Geometry.hh"

class Angular_Discretization;
class Energy_Discretization;
class Manufactured_Cross_Sections;
class Manufactured_Solution;
class Material_Factory;

/*
  Minimal class for using a solid geometry with a manufactured solution
*/
class Manufactured_Solid_Geometry : public Solid_Geometry
{
public:

    Manufactured_Solid_Geometry(std::shared_ptr<Angular_Discretization> angular,
                                std::shared_ptr<Energy_Discretization> energy,
                                std::shared_ptr<Manufactured_Solution> solution,
                                std::shared_ptr<Manufactured_Cross_Sections> cross_sections);

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
    virtual void output(XML_Node output_node) const override;

private:

    std::shared_ptr<Angular_Discretization> angular_;
    std::shared_ptr<Energy_Discretization> energy_;
    std::shared_ptr<Manufactured_Solution> solution_;
    std::shared_ptr<Manufactured_Cross_Sections> cross_sections_;
    std::shared_ptr<Material_Factory> material_factory_;
};

#endif

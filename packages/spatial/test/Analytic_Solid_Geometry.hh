#include <functional>

#include "Angular_Discretization.hh"
#include "Check.hh"
#include "Cross_Section.hh"
#include "Energy_Discretization.hh"
#include "Material.hh"
#include "Solid_Geometry.hh"
#include "Weak_Spatial_Discretization.hh"
#include "Weight_Function.hh"
#include "XML_Node.hh"

// All cross sections trivial except internal source, which is a given function
class Analytic_Solid_Geometry : public Solid_Geometry
{
public:
    Analytic_Solid_Geometry(int dimension,
                            std::shared_ptr<Angular_Discretization> angular,
                            std::shared_ptr<Energy_Discretization> energy,
                            std::function<double(std::vector<double> const &)> source):
        dimension_(dimension),
        angular_(angular),
        energy_(energy),
        source_(source)
    {
        Cross_Section::Dependencies none_group;
        none_group.energy = Cross_Section::Dependencies::Energy::GROUP;
        Cross_Section::Dependencies scattering_group2;
        scattering_group2.angular = Cross_Section::Dependencies::Angular::SCATTERING_MOMENTS;
        scattering_group2.energy = Cross_Section::Dependencies::Energy::GROUP_TO_GROUP;

        std::vector<double> one(1, 1.);
        std::vector<double> zero(1, 0.);
        sigma_t_ = std::make_shared<Cross_Section>(none_group,
                                                   angular_,
                                                   energy_,
                                                   one);
        sigma_s_ = std::make_shared<Cross_Section>(scattering_group2,
                                                   angular_,
                                                   energy_,
                                                   zero);
        nu_ = std::make_shared<Cross_Section>(none_group,
                                              angular_,
                                              energy_,
                                              one);
        sigma_f_ = std::make_shared<Cross_Section>(none_group,
                                                   angular_,
                                                   energy_,
                                                   zero);
        chi_ = std::make_shared<Cross_Section>(none_group,
                                               angular_,
                                               energy_,
                                               one);
    }

    virtual int dimension() const override
    {
        return dimension_;
    }
    virtual int find_region(std::vector<double> const &position) const override
    {
        AssertMsg(false, "not implemented");
        return -1;
    }
    virtual int find_surface(std::vector<double> const &position) const override
    {
        AssertMsg(false, "not implemented");
        return -1;
    }
    virtual int next_intersection(std::vector<double> const &initial_position,
                                  std::vector<double> const &initial_direction,
                                  int &final_region,
                                  double &distance,
                                  std::vector<double> &final_position) const override
    {
        AssertMsg(false, "not implemented");
        return -1;
    }
    virtual int next_boundary(std::vector<double> const &initial_position,
                              std::vector<double> const &initial_direction,
                              int &boundary_region,
                              double &distance,
                              std::vector<double> &final_position) const override
    {
        AssertMsg(false, "not implemented");
        return -1;
    }
    virtual void optical_distance(std::vector<double> const &initial_position,
                                  std::vector<double> const &final_position,
                                  std::vector<double> &optical_distance) const override
    {
        AssertMsg(false, "not implemented");
    }
    virtual std::shared_ptr<Material> material(std::vector<double> const &position) const override
    {
        Cross_Section::Dependencies none_group;
        none_group.energy = Cross_Section::Dependencies::Energy::GROUP;
        std::vector<double> data(1, source_(position));
        std::shared_ptr<Cross_Section> internal_source
            = std::make_shared<Cross_Section>(none_group,
                                              angular_,
                                              energy_,
                                              data);
        
        return std::make_shared<Material>(0,
                                          angular_,
                                          energy_,
                                          sigma_t_,
                                          sigma_s_,
                                          nu_,
                                          sigma_f_,
                                          chi_,
                                          internal_source);
        
    }
    virtual std::shared_ptr<Boundary_Source> boundary_source(int boundary_index,
                                                             std::vector<double> const &position) const override
    {
        AssertMsg(false, "not implemented");

        return std::shared_ptr<Boundary_Source>();
    }
    virtual void check_class_invariants() const override
    {
        AssertMsg(false, "not implemented");
    }
    virtual void output(XML_Node output_node) const override
    {
        AssertMsg(false, "not implemented");
    }

private:
    
    int dimension_;
    std::shared_ptr<Cross_Section> sigma_t_;
    std::shared_ptr<Cross_Section> sigma_s_;
    std::shared_ptr<Cross_Section> nu_;
    std::shared_ptr<Cross_Section> sigma_f_;
    std::shared_ptr<Cross_Section> chi_;
    std::shared_ptr<Angular_Discretization> angular_;
    std::shared_ptr<Energy_Discretization> energy_;
    std::function<double(std::vector<double> const &)> source_;
};

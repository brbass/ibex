#ifndef Weak_Spatial_Discretization_hh
#define Weak_Spatial_Discretization_hh

#include "Spatial_Discretization.hh"
#include "Weight_Function.hh"

class Basis_Function;
class KD_Tree;

class Weak_Spatial_Discretization : public Spatial_Discretization
{
public:

    // Constructor
    Weak_Spatial_Discretization(std::vector<std::shared_ptr<Basis_Function> > &bases,
                                std::vector<std::shared_ptr<Weight_Function> > &weights,
                                std::shared_ptr<Dimensional_Moments> dimensional_moments,
                                std::shared_ptr<KD_Tree> kd_tree = std::shared_ptr<KD_Tree>());

    // Point functions
    virtual bool has_reflection() const override
    {
        return has_reflection_;
    }
    virtual int number_of_points() const override
    {
        return number_of_points_;
    }
    virtual int number_of_boundary_points() const override
    {
        return number_of_boundary_bases_;
    }
    virtual int dimension() const override
    {
        return dimension_;
    }
    virtual int number_of_nodes() const override
    {
        return number_of_nodes_;
    }
    virtual int number_of_dimensional_moments() const override
    {
        return number_of_dimensional_moments_;
    }
    virtual std::shared_ptr<Dimensional_Moments> dimensional_moments() const override
    {
        return dimensional_moments_;
    }
    virtual std::shared_ptr<Point> point(int point_index) const override
    {
        return weights_[point_index];
    }
    virtual void output(XML_Node output_node) const override;
    virtual void check_class_invariants() const override;

    // Weight_Function functions
    virtual bool include_supg() const
    {
        return include_supg_;
    }
    virtual int number_of_boundary_weights() const
    {
        return number_of_boundary_weights_;
    }
    virtual std::vector<int> const &number_of_basis_functions() const
    {
        return number_of_basis_functions_;
    }
    virtual std::shared_ptr<Weight_Function> weight(int point_index) const
    {
        return weights_[point_index];
    }
    virtual std::shared_ptr<Basis_Function> basis(int point_index) const
    {
        return bases_[point_index];
    }
    virtual std::shared_ptr<Weight_Function> boundary_weight(int boundary_index) const
    {
        return boundary_weights_[boundary_index];
    }
    virtual std::shared_ptr<Basis_Function> boundary_basis(int boundary_index) const
    {
        return boundary_bases_[boundary_index];
    }
    virtual int nearest_point(std::vector<double> const &position) const;
    virtual double collocation_value(int i,
                                     std::vector<double> const &coefficients) const;
    virtual double weighted_collocation_value(int i,
                                              std::vector<double> const &coefficients) const;
    virtual void collocation_values(std::vector<double> const &coefficients,
                                    std::vector<double> &values) const;
    virtual void weighted_collocation_values(std::vector<double> const &coefficients,
                                             std::vector<double> &values) const;
    virtual double expansion_value(int i,
                                   std::vector<double> const &position,
                                   std::vector<double> const &coefficients) const;
    virtual double expansion_value(std::vector<double> const &position,
                                   std::vector<double> const &coefficients) const;
private:

    bool has_reflection_;
    bool include_supg_;
    int number_of_points_;
    int number_of_boundary_weights_;
    int number_of_boundary_bases_;
    int dimension_;
    int number_of_nodes_;
    int number_of_dimensional_moments_;
    std::vector<int> number_of_basis_functions_;
    std::vector<std::shared_ptr<Weight_Function> > weights_;
    std::vector<std::shared_ptr<Weight_Function> > boundary_weights_;
    std::vector<std::shared_ptr<Basis_Function> > bases_;
    std::vector<std::shared_ptr<Basis_Function> > boundary_bases_;
    std::shared_ptr<Dimensional_Moments> dimensional_moments_;
    std::shared_ptr<KD_Tree> kd_tree_;
    
};

#endif

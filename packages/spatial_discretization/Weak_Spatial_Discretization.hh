#ifndef Weak_Spatial_Discretization_hh
#define Weak_Spatial_Discretization_hh

#include "Spatial_Discretization.hh"
#include "Weight_Function.hh"

class Basis_Function;

class Weak_Spatial_Discretization : public Spatial_Discretization
{
public:

    // Constructor
    Weak_Spatial_Discretization(std::vector<std::shared_ptr<Basis_Function> > &bases,
                                std::vector<std::shared_ptr<Weight_Function> > &weights);

    virtual int number_of_points() const override
    {
        return number_of_points_;
    }
    virtual int number_of_boundary_points() const override
    {
        return number_of_boundary_weights_;
    }
    virtual int dimension() const override
    {
        return dimension_;
    }
    virtual int number_of_nodes() const override
    {
        return number_of_nodes_;
    }
    virtual std::shared_ptr<Point> point(int point_index) const override
    {
        return weights_[point_index];
    }
    virtual std::shared_ptr<Weight_Function> weight(int point_index) const
    {
        return weights_[point_index];
    }
    // virtual std::shared_ptr<Basis_Function> basis(int point_index) const
    // {
    //     return bases_[point_index];
    // }
    virtual void output(XML_Node output_node) const override;
    virtual void check_class_invariants() const override;

private:

    int number_of_points_;
    int number_of_boundary_weights_;
    int number_of_boundary_bases_;
    int dimension_;
    int number_of_nodes_;
    std::vector<std::shared_ptr<Weight_Function> > weights_;
    std::vector<std::shared_ptr<Weight_Function> > boundary_weights_;
    std::vector<std::shared_ptr<Basis_Function> > bases_;
    std::vector<std::shared_ptr<Basis_Function> > boundary_bases_;
};

#endif

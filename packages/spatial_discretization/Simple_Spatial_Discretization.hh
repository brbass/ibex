#ifndef Simple_Spatial_Discretization_hh
#define Simple_Spatial_Discretization_hh

#include "Spatial_Discretization.hh"

#include "Point.hh"

class Simple_Spatial_Discretization : public Spatial_Discretization
{
public:
    
    Simple_Spatial_Discretization(vector<shared_ptr<Point> > points);
    
    virtual int number_of_points() const override
    {
        return number_of_points_;
    }

    virtual int number_of_internal_points() const override
    {
        return number_of_internal_points_;
    }

    virtual int number_of_boundary_points() const override
    {
        return number_of_boundary_points_;
    }
    
    virtual int dimension() const override
    {
        return dimension_;
    }

    virtual int number_of_nodes() const override
    {
        return 1;
    }

    virtual vector<int> const &boundary_points() const override
    {
        return boundary_points_;
    }
    
    virtual vector<int> const &internal_points() const override
    {
        return internal_points_;
    }

    virtual shared_ptr<Point> point(int point_index) const override
    {
        return points_[point_index];
    }
    
    virtual void output(pugi::xml_node &output_node) const override;

    virtual void check_class_invariants() const override;
    
private:

    int dimension_;
    int number_of_points_;
    int number_of_internal_points_;
    int number_of_boundary_points_;
    vector<int> boundary_points_;
    vector<int> internal_points_;
    vector<shared_ptr<Point> > points_;
};

#endif

#ifndef Constructive_Solid_Geometry_hh
#define Constructive_Solid_Geometry_hh

#include "Boundary_Source.hh"
#include "Material.hh"
#include "Region.hh"
#include "Solid_Geometry.hh"
#include "Surface.hh"

/* 
   Solid geometry for Monte Carlo and RBF calculations
   All methods should allow for final and initial values of the same type
   to be the same vectors
*/
class Constructive_Solid_Geometry : public Solid_Geometry
{
public:

    Constructive_Solid_Geometry(int dimension,
                                vector<shared_ptr<Surface> > const &surfaces,
                                vector<shared_ptr<Region> > const &regions,
                                vector<shared_ptr<Material> > const &materials,
                                vector<shared_ptr<Boundary_Source> > const &boundary_sources);
    
    virtual bool cartesian_boundaries() const
    {
        return cartesian_boundaries_;
    }
    virtual int dimension() const override
    {
        return dimension_;
    }
    virtual int number_of_surfaces() const
    {
        return surfaces_.size();
    }
    virtual int number_of_regions() const
    {
        return regions_.size();
    }
    virtual double delta_distance() const
    {
        return delta_distance_;
    }
    
    virtual shared_ptr<Surface> const &surface(int s) const
    {
        return surfaces_[s];
    }
    virtual shared_ptr<Region> const &region(int r) const
    {
        return regions_[r];
    }
    
    virtual int find_region(vector<double> const &position) const override;
    virtual int find_region_including_surface(vector<double> const &position) const;
    virtual int find_surface(vector<double> const &position) const override;

    virtual int next_intersection(vector<double> const &initial_position,
                                  vector<double> const &initial_direction,
                                  int &final_region,
                                  double &distance,
                                  vector<double> &final_position) const override;
    virtual int next_intersection(int initial_region,
                                  vector<double> const &initial_position,
                                  vector<double> const &initial_direction,
                                  int &final_region,
                                  double &distance,
                                  vector<double> &final_position) const;
    
    virtual int next_boundary(vector<double> const &initial_position,
                              vector<double> const &initial_direction,
                              int &boundary_region,
                              double &distance,
                              vector<double> &final_position) const override;
    virtual int next_boundary(int initial_region,
                              vector<double> const &initial_position,
                              vector<double> const &initial_direction,
                              int &boundary_region,
                              double &distance,
                              vector<double> &final_position) const;
    
    virtual void new_position(double distance,
                              vector<double> const &initial_position,
                              vector<double> const &initial_direction,
                              vector<double> &final_position) const;
    
    virtual void optical_distance(vector<double> const &initial_position,
                                  vector<double> const &final_position,
                                  vector<double> &optical_distance) const override;
    virtual shared_ptr<Material> material(vector<double> const &position) const override;
    
    virtual void check_class_invariants() const override;
    virtual void output(pugi::xml_node &output_node) const override;
    
protected:

    bool cartesian_boundaries_;
    int dimension_;
    double optical_tolerance_;
    double delta_distance_;
    vector<shared_ptr<Surface> > boundary_surfaces_;
    vector<shared_ptr<Surface> > surfaces_;
    vector<shared_ptr<Region> > regions_;
    vector<shared_ptr<Material> > materials_;
    vector<shared_ptr<Boundary_Source> > boundary_sources_;
};

#endif

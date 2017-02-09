#ifndef Constructive_Solid_Geometry_hh
#define Constructive_Solid_Geometry_hh

#include "Solid_Geometry.hh"

class Boundary_Source;
class Region;
class Surface;

/* 
   Solid geometry for Monte Carlo and RBF calculations
   All methods should allow for final and initial values of the same type
   to be the same vectors
*/
class Constructive_Solid_Geometry : public Solid_Geometry
{
public:

    Constructive_Solid_Geometry(int dimension,
                                std::vector<std::shared_ptr<Surface> > const &surfaces,
                                std::vector<std::shared_ptr<Region> > const &regions,
                                std::vector<std::shared_ptr<Material> > const &materials,
                                std::vector<std::shared_ptr<Boundary_Source> > const &boundary_sources);
    
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
    virtual int number_of_boundary_surfaces() const
    {
        return boundary_surfaces_.size();
    }
    virtual double delta_distance() const
    {
        return delta_distance_;
    }
    
    virtual std::shared_ptr<Surface> surface(int s) const
    {
        return surfaces_[s];
    }
    virtual std::shared_ptr<Surface> boundary_surface(int s) const
    {
        return boundary_surfaces_[s];
    }
    virtual std::shared_ptr<Cartesian_Plane> cartesian_boundary_surface(int s) const
    {
        return cartesian_boundary_surfaces_[s];
    }
    virtual std::vector<std::shared_ptr<Surface> > surfaces() const
    {
        return surfaces_;
    }
    virtual std::vector<std::shared_ptr<Surface> > boundary_surfaces() const
    {
        return boundary_surfaces_;
    }
    virtual std::vector<std::shared_ptr<Cartesian_Plane> > boundary_surfaces() const
    {
        return cartesian_boundary_surfaces_;
    }
    virtual std::shared_ptr<Region> region(int r) const
    {
        return regions_[r];
    }
    
    virtual int find_region(std::vector<double> const &position) const override;
    virtual int find_region_including_surface(std::vector<double> const &position) const;
    virtual int find_surface(std::vector<double> const &position) const override;

    virtual int next_intersection(std::vector<double> const &initial_position,
                                  std::vector<double> const &initial_direction,
                                  int &final_region,
                                  double &distance,
                                  std::vector<double> &final_position) const override;
    virtual int next_intersection(int initial_region,
                                  std::vector<double> const &initial_position,
                                  std::vector<double> const &initial_direction,
                                  int &final_region,
                                  double &distance,
                                  std::vector<double> &final_position) const;
    
    virtual int next_boundary(std::vector<double> const &initial_position,
                              std::vector<double> const &initial_direction,
                              int &boundary_region,
                              double &distance,
                              std::vector<double> &final_position) const override;
    virtual int next_boundary(int initial_region,
                              std::vector<double> const &initial_position,
                              std::vector<double> const &initial_direction,
                              int &boundary_region,
                              double &distance,
                              std::vector<double> &final_position) const;
    
    virtual void new_position(double distance,
                              std::vector<double> const &initial_position,
                              std::vector<double> const &initial_direction,
                              std::vector<double> &final_position) const;
    
    virtual void optical_distance(std::vector<double> const &initial_position,
                                  std::vector<double> const &final_position,
                                  std::vector<double> &optical_distance) const override;
    virtual std::shared_ptr<Material> material(std::vector<double> const &position) const override;
    
    virtual void check_class_invariants() const override;
    virtual void output(XML_Node output_node) const override;
    
protected:

    bool cartesian_boundaries_;
    int dimension_;
    double optical_tolerance_;
    double delta_distance_;
    std::vector<std::shared_ptr<Surface> > surfaces_;
    std::vector<std::shared_ptr<Surface> > boundary_surfaces_;
    std::vector<std::shared_ptr<Cartesian_Plane> > cartesian_boundary_surfaces_;
    std::vector<std::shared_ptr<Region> > regions_;
    std::vector<std::shared_ptr<Material> > materials_;
    std::vector<std::shared_ptr<Boundary_Source> > boundary_sources_;
};

#endif

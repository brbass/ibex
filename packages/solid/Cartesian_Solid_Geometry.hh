#ifndef Cartesian_Solid_Geometry_hh
#define Cartesian_Solid_Geometry_hh

#include "Solid_Geometry.hh"

/* 
   Cartesian solid geometry with equal cell length for all cells in each dimension
   Allows easy cell finding and boundary checking
   Material index is i_x + n_x * (i_y + n_y * i_z)
*/
class Cartesian_Solid_Geometry : public Solid_Geometry
{
public:
    
    Cartesian_Solid_Geometry(int dimension,
                             std::vector<int> const &number_of_cells,
                             std::vector<double> const &total_length,
                             std::vector<std::shared_ptr<Material> > const &materials);
    
    virtual int dimension() const override
    {
        return dimension_;
    }
    virtual int find_region(std::vector<double> const &position) const override;
    virtual int find_surface(std::vector<double> const &position) const override;
    virtual int next_intersection(std::vector<double> const &initial_position,
                                  std::vector<double> const &initial_direction,
                                  int &final_region,
                                  double &distance,
                                  std::vector<double> &final_position) const override;
    virtual int next_boundary(std::vector<double> const &initial_position,
                              std::vector<double> const &initial_direction,
                              int &boundary_region,
                              double &distance,
                              std::vector<double> &final_position) const override;
    virtual void optical_distance(std::vector<double> const &initial_position,
                                  std::vector<double> const &final_position,
                                  std::vector<double> &optical_distance) const override;
    virtual std::shared_ptr<Material> material(std::vector<double> const &position) const override;
    virtual std::shared_ptr<Boundary_Source> boundary_source(std::vector<double> const &position) const override;
    virtual void check_class_invariants() const override;
    virtual void output(XML_Node output_node) const override;
    
private:

    int dimension_;
    std::vector<int> number_of_cells_;
    std::vector<double> total_length_;
    std::vector<double> cell_length_;
    std::vector<std::shared_ptr<Material> > materials_;
};

#endif

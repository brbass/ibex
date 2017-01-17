#ifndef Cartesian_Overlay_hh
#define Cartesian_Overlay_hh

#include <memory>
#include <vector>

#include "Indexing.hh"

class Distance;

class Cartesian_Overlay
{
public:
    
    enum Cell_Error
    {
        LOWER_BOUND = -1,
        UPPER_BOUND = -2
    };
    
    Cartesian_Overlay(int dimension,
                      std::vector<int> const &number_of_cells_per_dimension,
                      std::vector<double> const &cartesian_bounds,
                      std::shared_ptr<Distance> distance);
    
    int number_of_points() const
    {
        return number_of_points_;
    }
    void add_point(std::vector<double> const &position);
    int get_cell(std::vector<double> const &position) const;
    // void nearest_points(int number_of_points,
    //                     std::vector<double> const &position,
    //                     std::vector<int> &indices,
    //                     std::vector<double> &distances) const;
    bool has_neighbor(double minimum_distance,
                      std::vector<double> const &position) const;
    void get_cells_to_check(int cell,
                            std::vector<int> const &number_of_cells_to_check,
                            std::vector<int> &cells_to_check) const;

    virtual void check_class_invariants() const;
    
private:
    
    // Static variables
    int dimension_;
    int number_of_cells_;
    std::vector<int> number_of_cells_per_dimension_;
    std::vector<double> cell_length_;
    std::vector<double> cartesian_bounds_;
    std::shared_ptr<Distance> distance_;
    std::shared_ptr<Indexing<int> > indexing_;
    
    // Dynamic variables
    int number_of_points_;
    std::vector<int> cell_by_index_; // which cell each point is in
    std::vector<std::vector<int> > indices_by_cell_; // points that are in each cell
    std::vector<std::vector<double> > position_by_index_;
};

#endif

#ifndef Cartesian_Overlay_hh
#define Cartesian_Overlay_hh

#include <memory>
#include <vector>

#include "Indexing.hh"

using std::shared_ptr;
using std::vector;

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
                      vector<int> const &number_of_cells_per_dimension,
                      vector<double> const &cartesian_bounds,
                      shared_ptr<Distance> distance);
    
    int number_of_points() const
    {
        return number_of_points_;
    }
    void add_point(vector<double> const &position);
    int get_cell(vector<double> const &position) const;
    // void nearest_points(int number_of_points,
    //                     vector<double> const &position,
    //                     vector<int> &indices,
    //                     vector<double> &distances) const;
    bool has_neighbor(double minimum_distance,
                      vector<double> const &position) const;
    void get_cells_to_check(int cell,
                            vector<int> const &number_of_cells_to_check,
                            vector<int> &cells_to_check) const;

    virtual void check_class_invariants() const;
    
private:
    
    // Static variables
    int dimension_;
    int number_of_cells_;
    vector<int> number_of_cells_per_dimension_;
    vector<double> cell_length_;
    vector<double> cartesian_bounds_;
    shared_ptr<Distance> distance_;
    shared_ptr<Indexing<int> > indexing_;
    
    // Dynamic variables
    int number_of_points_;
    vector<int> cell_by_index_; // which cell each point is in
    vector<vector<int> > indices_by_cell_; // points that are in each cell
    vector<vector<double> > position_by_index_;
};

#endif

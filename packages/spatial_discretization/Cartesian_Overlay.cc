#include "Cartesian_Overlay.hh"

#include <cmath>

#include "Check.hh"
#include "Distance.hh"
#include "Indexing.hh"

using std::floor;
using std::make_shared;

Cartesian_Overlay::
Cartesian_Overlay(int dimension,
                  vector<int> const &number_of_cells_per_dimension,
                  vector<double> const &cartesian_bounds,
                  shared_ptr<Distance> distance):
    dimension_(dimension),
    number_of_cells_per_dimension_(number_of_cells_per_dimension),
    cartesian_bounds_(cartesian_bounds),
    distance_(distance)
{
    indexing_ = make_shared<Indexing<int> >(dimension,
                                            number_of_cells_per_dimension);
        
    cell_length_.resize(dimension);
    
    for (int d = 0; d < dimension_; ++d)
    {
        cell_length_[d] = (cartesian_bounds[d + dimension * 1] - cartesian_bounds[d + dimension * 0]) / number_of_cells_per_dimension[d];
    }

    number_of_cells_ = 1;

    for (int d = 0; d < dimension_; ++d)
    {
        number_of_cells_ *= number_of_cells_per_dimension_[d];
    }
    
    number_of_points_ = 0;
    cell_by_index_.resize(0);
    indices_by_cell_.assign(number_of_cells_, vector<int>());
    position_by_index_.resize(0);

    check_class_invariants();
}

void Cartesian_Overlay::
add_point(vector<double> const &position)
{
    int index = number_of_points_;
    int cell = get_cell(position);
    
    number_of_points_ += 1;
    position_by_index_.push_back(position);
    cell_by_index_.push_back(cell);
    indices_by_cell_[cell].push_back(index);
}

int Cartesian_Overlay::
get_cell(vector<double> const &position) const
{
    vector<int> subscript(dimension_);
    
    for (int d = 0; d < dimension_; ++d)
    {
        double x = position[d] - cartesian_bounds_[d + dimension_ * 0];
        double y = cartesian_bounds_[d + dimension_ * 1] - position[d];

        int c = floor(x / cell_length_[d]);
        
        if (c < 0)
        {
            subscript[d] = LOWER_BOUND;
            
            AssertMsg(false, "point outside Cartesian bounds");
        }
        else if (c > number_of_cells_per_dimension_[d])
        {
            subscript[d] = UPPER_BOUND;
            
            AssertMsg(false, "point outside Cartesian bounds");
        }
        else
        {
            subscript[d] = c;
        }
    }
    
    return indexing_->subscript_to_index(subscript);
}

bool Cartesian_Overlay::
has_neighbor(double minimum_distance,
             vector<double> const &position) const
{
    int cell = get_cell(position);
    vector<int> number_of_cells_to_check(dimension_);
    for (int d = 0; d < dimension_; ++d)
    {
        number_of_cells_to_check[d] = ceil(minimum_distance / cell_length_[d]);
    }
    
    vector<int> cells_to_check;
    get_cells_to_check(cell,
                       number_of_cells_to_check,
                       cells_to_check);
    
    int total_number_of_cells = cells_to_check.size();
    
    for (int i = 0; i < total_number_of_cells; ++i)
    {
        int j = cells_to_check[i]; // spatial cell
        vector<int> const indices = indices_by_cell_[j]; // point indices
        int number_of_points_in_cell = indices.size();

        for (int k = 0; k < number_of_points_in_cell; ++k)
        {
            vector<double> const cell_position = position_by_index_[indices[k]];
            
            double distance = distance_->distance(position,
                                                  cell_position);

            if (distance < minimum_distance)
            {
                return true;
            }
        }
    }
    
    return false;
}

void Cartesian_Overlay::
get_cells_to_check(int cell,
                   vector<int> const &number_of_cells_to_check,
                   vector<int> &cells_to_check) const
{
    cells_to_check.resize(0);
    vector<int> subscript = indexing_->index_to_subscript(cell);
    vector<int> starting_subscript(dimension_);
    vector<int> ending_subscript(dimension_);
    for (int d = 0; d < dimension_; ++d)
    {
        starting_subscript[d] = subscript[d] - number_of_cells_to_check[d];
        ending_subscript[d] = subscript[d] + number_of_cells_to_check[d];
        
        if (starting_subscript[d] < 0)
        {
            starting_subscript[d] = 0;
        }
        if (ending_subscript[d] > number_of_cells_per_dimension_[d] - 1)
        {
            ending_subscript[d] = number_of_cells_per_dimension_[d] - 1;
        }
    }
    
    bool completed = false;
    vector<int> current_subscript(starting_subscript);
    while (!completed)
    {
        cells_to_check.push_back(indexing_->subscript_to_index(current_subscript));
        
        for (int d = 0; d < dimension_; ++d)
        {
            current_subscript[d] += 1;
            
            if (current_subscript[d] <= ending_subscript[d])
            {
                break;
            }
            else
            {
                if (d == dimension_ - 1)
                {
                    completed = true;
                    break;
                }
                else
                {
                    current_subscript[d] = starting_subscript[d];
                }
            }
        }
    }
}

void Cartesian_Overlay::
check_class_invariants() const
{
    Assert(number_of_cells_per_dimension_.size() == dimension_);
    Assert(cell_length_.size() == dimension_);
    Assert(cartesian_bounds_.size() == dimension_ * 2);
    Assert(distance_);
    Assert(indexing_);

    Assert(cell_by_index_.size() == number_of_points_);
    Assert(indices_by_cell_.size() == number_of_cells_);
    Assert(position_by_index_.size() == number_of_points_);
}

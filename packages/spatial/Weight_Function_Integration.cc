#include "Weight_Function_Integration.hh"

using std::shared_ptr;
using std::vector;

Weight_Function::
Weight_Function(Weight_Function::Options options,
                vector<shared_ptr<Meshless_Function> > const &bases,
                vector<shared_ptr<Meshless_Function> > const &weights,
                shared_ptr<Solid_Geometry> solid,
                vector<vector<double> > limits,
                vector<int> num_intervals):
    options_(options),
    bases_(bases),
    weights_(weights),
    solid_(solid),
    mesh_(solid->dimension(),
          limits,
          num_intervals)
{
}

Weight_Function::Mesh
Mesh(int dimension,
     vector<vector<double> > limits,
     vector<int> num_intervals):
    dimension_(dimension),
    limits_(limits),
    num_intervals_(num_intervals)
{
    kd_tree_ = get_kd_tree();
    
    
}

shared_ptr<KD_Tree> Weight_Function::Mesh::
get_kd_tree() const
{
    vector<vector<double> > points;
    
    // Check sizes
    Assert(num_intervals_.size() == dimension);
    Assert(limits_.size() == dimension);
    
    // Get total number of points
    number_of_points_ = 1;
    number_of_cells_ = 1;
    for (int d = 0; d < dimension_; ++d)
    {
        Assert(num_intervals_[d] >= 1);
        number_of_points_ *= num_intervals_[d] + 1;
        number_of_cells_ *= num_intervals_[d];
    }
    
    // Get Cartesian grid of points
    points.assign(number_of_points_, vector<double>(dimension_));
    int index = 0;
    vector<int> indices(dimension, 0); // current index
    intervals_.resize(dimension);
    for (int d = 0; d < dimension; ++d)
    {
        intervals_[d] = (limits_[d][1] - limits_[d][0]) / static_cast<double>(num_intervals_[d]);
    }
    while (index < number_of_points_)
    {
        for (int d = 0; d < dimension; ++d)
        {
            points[index][d] = limits_[d][0] + intervals_[d] * indices[d];
        }
        
        index += 1;
        indices[0] += 1;
        for (int d = 0; d < dimension - 1; ++d)
        {
            if (indices[d] > num_intervals_[d])
            {
                indices[d] = 0;
                indices[d + 1] += 1;
            }
        }
    }

    return make_shared<KD_Tree>(dimension_,
                                number_of_points_,
                                points);
}

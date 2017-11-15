#include "KD_Tree.hh"

#include "Check.hh"

using namespace std;

KD_Tree::
KD_Tree(int dimension,
        int number_of_points,
        vector<vector<double> > const &points):
    dimension_(dimension),
    number_of_points_(number_of_points),
    points_(points)
{
    Assert(dimension >= 1);
    
    adaptor_ = make_shared<KD_Adaptor>(*this);
    kd_tree_ = make_shared<KDT>(dimension_,
                                *adaptor_);
    kd_tree_->buildIndex();
}

KD_Tree::KD_Adaptor::
KD_Adaptor(KD_Tree const &tree):
    tree_(tree)
{
}

void KD_Tree::
find_neighbors(int number_of_neighbors,
               vector<double> const &position,
               vector<int> &indices,
               vector<double> &squared_distances) const
{
    Check(number_of_neighbors <= number_of_points_);
    
    indices.resize(number_of_neighbors);
    squared_distances.resize(number_of_neighbors);

    kd_tree_->knnSearch(&position[0],
                        number_of_neighbors,
                        &indices[0],
                        &squared_distances[0]);
}

int KD_Tree::
radius_search(double radius,
              vector<double> const &position,
              vector<int> &indices,
              vector<double> &squared_distances) const
{
    vector<pair<int, double> > indices_and_distances;
    nanoflann::SearchParams search_parameters;
    search_parameters.sorted = true;
    
    int number_of_neighbors = kd_tree_->radiusSearch(&position[0],
                                                     radius*radius,
                                                     indices_and_distances,
                                                     search_parameters);

    indices.resize(number_of_neighbors);
    squared_distances.resize(number_of_neighbors);
    for (int i = 0; i < number_of_neighbors; ++i)
    {
        indices[i] = indices_and_distances[i].first;
        squared_distances[i] = indices_and_distances[i].second;
    }
    
    return number_of_neighbors;
}

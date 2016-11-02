#include "KD_Tree.hh"

#include "Check.hh"

using namespace std;

KD_Tree::
KD_Tree(int dimension,
        int number_of_points,
        vector<double> const &points):
    dimension_(dimension),
    number_of_points_(number_of_points),
    points_(points)
{
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
find_neighbors(int index,
               int number_of_neighbors,
               vector<int> &indices,
               vector<double> &distances) const
{
    Check(number_of_neighbors <= number_of_points_);
    
    indices.resize(number_of_neighbors);
    distances.resize(number_of_neighbors);

    vector<double> points(dimension_);
    
    kd_tree_->knnSearch(&points_[dimension_ * index],
                        number_of_neighbors,
                        &indices[0],
                        &distances[0]);
}


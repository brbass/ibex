#ifndef KD_Tree_hh
#define KD_Tree_hh

#include <memory>
#include <vector>

#include "nanoflann.hh"

/*
  KD Tree for finding nearest neighbors
  General adaptor for nanoflann
*/
class KD_Tree
{
    class KD_Adaptor;
    
    // Distance adaptor
    typedef nanoflann::L2_Adaptor<double,
                                  KD_Adaptor> L2A;
    // KD tree type
    typedef nanoflann::KDTreeSingleIndexAdaptor<L2A,
                                                KD_Adaptor,
                                                -1,
                                                int> KDT;
    
public:

    // Constructor
    KD_Tree(int dimension,
            int number_of_points,
            std::vector<std::vector<double> > const &points);

    // Find nearest neighbor indices and distances given a position
    virtual void find_neighbors(int number_of_neighbors,
                                std::vector<double> const &position,
                                std::vector<int> &indices,
                                std::vector<double> &distances) const;
    
    // Find all points within a radius of the position; return number of matches
    virtual int radius_search(double radius,
                              std::vector<double> const &position,
                              std::vector<int> &indices,
                              std::vector<double> &distances) const;

    // Get point
    virtual std::vector<double> const &point(int i) const
    {
        return points_[i];
    }
    
private:
    
    int number_of_points_;
    int dimension_;
    
    std::vector<std::vector<double> > points_;
    
    std::shared_ptr<KD_Adaptor> adaptor_;
    std::shared_ptr<KDT> kd_tree_;
    
    class KD_Adaptor
    {
    public:

        KD_Adaptor(KD_Tree const &tree);
        
        inline int kdtree_get_point_count() const
        {
            return tree_.number_of_points_;
        }
        
        inline double kdtree_get_pt(const int idx, int dim) const
        {
            return tree_.points_[idx][dim];
        }
        
        inline double kdtree_distance(const double *p1, const int idx_p2, int /*size*/) const
        {
            double sum = 0;
            
            for (int d = 0; d < tree_.dimension_; ++d)
            {
                double dist = (p1[d] - tree_.points_[idx_p2][d]);
                
                sum += dist * dist;
            }
            
            return sum;
        }
        
        template <class BBOX>
        bool kdtree_get_bbox(BBOX &/*bb*/) const
        {
            return false;
        }

    private:

        KD_Tree const &tree_;
    };
};

#endif

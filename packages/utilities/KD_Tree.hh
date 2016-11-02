#ifndef KD_Tree_hh
#define KD_Tree_hh

#include <memory>
#include <vector>

#include "nanoflann.hh"

class KD_Tree
{
    class KD_Adaptor;
    typedef nanoflann::L2_Adaptor<double,
                                  KD_Adaptor> L2A;
    typedef nanoflann::KDTreeSingleIndexAdaptor<L2A,
                                                KD_Adaptor,
                                                -1,
                                                int> KDT;
    
public:
    
    KD_Tree(int dimension,
            int number_of_points,
            std::vector<double> const &points);
    
    virtual void find_neighbors(int index,
                                int number_of_neighbors,
                                std::vector<int> &indices,
                                std::vector<double> &distances) const;
    
private:
    
    int number_of_points_;
    int dimension_;
    
    std::vector<double> points_;
    
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
            return tree_.points_[dim + tree_.dimension_ * idx];
        }
        
        inline double kdtree_distance(const double *p1, const int idx_p2, int /*size*/) const
        {
            double sum = 0;
            
            for (int i = 0; i < tree_.dimension_; ++i)
            {
                double dist = (p1[i] - tree_.points_[i + tree_.dimension_ * idx_p2]);
                
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

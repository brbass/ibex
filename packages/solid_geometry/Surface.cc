#include "Surface.hh"

#include <limits>

#include "Boundary_Source.hh"
#include "Vector_Functions_2D.hh"
#include "Vector_Functions_3D.hh"

using std::numeric_limits;

namespace vf2 = Vector_Functions_2D;
namespace vf3 = Vector_Functions_3D;

Surface::
Surface(int index,
        int dimension,
        Surface_Type surface_type):
    index_(index),
    dimension_(dimension),
    surface_type_(surface_type)
{
    intersection_tolerance_ = 10 * numeric_limits<double>::epsilon();
    relation_tolerance_ = 1000 * numeric_limits<double>::epsilon();
    normal_tolerance_ = 1000 * numeric_limits<double>::epsilon();
};

bool Surface::
reflected_direction(vector<double> const &position,
                    vector<double> const &old_direction,
                    vector<double> &new_direction,
                    bool check_normal)
{
    vector<double> normal;

    if(normal_direction(position,
                        normal,
                        check_normal))
    {
        switch(dimension_)
        {
        case 1:
            new_direction.assign(1, -old_direction[0]);
            
            break;
        case 2:
            new_direction = vf2::subtract(old_direction,
                                          vf2::multiply(normal,
                                                        2 * vf2::dot(old_direction,
                                                                     normal)));
            
            break;
        case 3:
            new_direction = vf3::subtract(old_direction,
                                          vf3::multiply(normal,
                                                        2 * vf3::dot(old_direction,
                                                                     normal)));
            
            break;
        }
        
        return true;
    }
    else
    {
        return false;
    }
}


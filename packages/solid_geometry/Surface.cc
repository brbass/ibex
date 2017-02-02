#include "Surface.hh"

#include <limits>

#include "Boundary_Source.hh"
#include "Vector_Functions_2D.hh"
#include "Vector_Functions_3D.hh"

using namespace std;

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

Surface::Reflection Surface::
reflected_direction(vector<double> const &position,
                    vector<double> const &old_direction,
                    vector<double> &new_direction,
                    bool check_normal)
{
    Normal normal = normal_direction(position,
                                     check_normal);
    Surface::Reflection reflection;
    if(normal.exists)
    {
        reflection.exists = true;
        switch(dimension_)
        {
        case 1:
            reflection.direction.assign(1, -old_direction[0]);
            
            break;
        case 2:
            reflection.direction = vf2::subtract(old_direction,
                                                 vf2::multiply(normal,
                                                               2 * vf2::dot(old_direction,
                                                                            normal)));
            
            break;
        case 3:
            reflection.direction = vf3::subtract(old_direction,
                                                 vf3::multiply(normal,
                                                               2 * vf3::dot(old_direction,
                                                                            normal)));
            
            break;
        }
        return reflection;
    }
    else
    {
        reflection.exists = false;
        return reflection;
    }
}


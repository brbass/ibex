#include "Point.hh"

Point::
Point()
{
}

string Point::
point_type_string() const
{
    switch(point_type())
    {
    case Point_Type::INTERNAL:
        return "internal";
    case Point_Type::BOUNDARY:
        return "boundary";
    }
}


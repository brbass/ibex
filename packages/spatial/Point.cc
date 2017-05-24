#include "Point.hh"

#include "Conversion.hh"

using std::make_shared;
using std::pair;
using std::shared_ptr;
using std::string;
using std::vector;

Point::
Point()
{
}

shared_ptr<Conversion<Point::Point_Type, string> > Point::
point_type_conversion() const
{
    vector<pair<Point_Type, string> > conversions
        = {{Point_Type::BOUNDARY, "boundary"},
           {Point_Type::INTERNAL, "internal"}};
    return make_shared<Conversion<Point_Type, string> >(conversions);
}

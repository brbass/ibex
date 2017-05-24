#include "Weighting_Operator.hh"

#include "Angular_Discretization.hh"
#include "Conversion.hh"
#include "Energy_Discretization.hh"
#include "Weak_Spatial_Discretization.hh"

using std::make_shared;
using std::pair;
using std::shared_ptr;
using std::string;
using std::vector;

Weighting_Operator::
Weighting_Operator(shared_ptr<Weak_Spatial_Discretization> spatial,
                   shared_ptr<Angular_Discretization> angular,
                   shared_ptr<Energy_Discretization> energy,
                   Options options):
    spatial_(spatial),
    angular_(angular),
    energy_(energy),
    options_(options)
{
}

std::shared_ptr<Conversion<Weighting_Operator::Options::Normalization, string> > Weighting_Operator::Options::
normalization_conversion() const
{
    vector<pair<Normalization, string> > conversions
        = {{Normalization::AUTO, "auto"},
           {Normalization::TRUE, "true"},
           {Normalization::FALSE, "false"}};
    return make_shared<Conversion<Normalization, string> >(conversions);
}

std::shared_ptr<Conversion<Weighting_Operator::Options::Include_SUPG, string> > Weighting_Operator::Options::
include_supg_conversion() const
{
    vector<pair<Include_SUPG, string> > conversions
        = {{Include_SUPG::AUTO, "auto"},
           {Include_SUPG::TRUE, "true"},
           {Include_SUPG::FALSE, "false"}};
    return make_shared<Conversion<Include_SUPG, string> >(conversions);
}

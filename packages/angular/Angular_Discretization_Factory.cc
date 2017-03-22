#include "Angular_Discretization_Factory.hh"

#include "Angular_Discretization.hh"
#include "Gauss_Legendre_Quadrature.hh"
#include "LDFE_Quadrature.hh"

using namespace std;

Angular_Discretization_Factory::
Angular_Discretization_Factory()
{
}

shared_ptr<Angular_Discretization> Angular_Discretization_Factory::
get_angular_discretization(int dimension,
                           int number_of_moments,
                           int angular_rule)
{
    if (dimension == 1)
    {
        return make_shared<Gauss_Legendre_Quadrature>(dimension,
                                                      number_of_moments,
                                                      angular_rule);
    }
    else
    {
        return make_shared<LDFE_Quadrature>(dimension,
                                            number_of_moments,
                                            angular_rule);
    }
}

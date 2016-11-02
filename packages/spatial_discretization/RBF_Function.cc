#include "RBF_Function.hh"

#include "Distance.hh"
#include "RBF.hh"
#include "XML_Functions.hh"

RBF_Function::
RBF_Function(shared_ptr<RBF> rbf,
             shared_ptr<Distance> distance):
    rbf_(rbf),
    distance_(distance)
{
}

double RBF_Function::
basis(double shape,
      vector<double> const &r,
      vector<double> const &r0) const
{
    double dist = distance_->distance(r,
                                      r0);
    return rbf_->basis(shape * dist);
}

double RBF_Function::
d_basis(int dim,
        double shape,
        vector<double> const &r,
        vector<double> const &r0) const
{
    double dist = distance_->distance(r,
                                      r0);
    double d_dist = distance_->d_distance(dim,
                                          r,
                                          r0);
    
    return rbf_->d_basis(shape * dist,
                         shape * d_dist);
}


double RBF_Function::
dd_basis(int dim,
         double shape,
         vector<double> const &r,
         vector<double> const &r0) const
{
    double dist = distance_->distance(r,
                                      r0);
    double d_dist = distance_->d_distance(dim,
                                          r,
                                          r0);
    double dd_dist = distance_->dd_distance(dim,
                                            r,
                                            r0);

    return rbf_->dd_basis(shape * dist,
                          shape * shape * d_dist * d_dist,
                          shape * dd_dist);
}

vector<double> RBF_Function::
gradient_basis(double shape,
               vector<double> const &r,
               vector<double> const &r0) const
{
    double dist = distance_->distance(group,
                                      r,
                                      r0);
    vector<double> grad = distance_->gradient_distance(r,
                                                       r0);

    int dimension = distance_->dimension();
    
    vector<double> res(dimension);
    
    for (int d = 0; d < dimension; ++d)
    {
        res[d] = rbf_->d_basis(shape * dist,
                               shape * grad[d]);
    }

    return res;
}

double RBF_Function::
laplacian(double shape,
          vector<double> const &r,
          vector<double> const &r0) const
{
    double dist = distance_->distance(r,
                                      r0);
    vector<double> grad = distance_->gradient_distance(r,
                                                       r0);
    double lap = distance_->laplacian_distance(r,
                                               r0);

    int dimension = distance_->dimension();

    double grad2 = 0;
    for (int d = 0; d < dimension; ++d)
    {
        grad2 += grad[d] * grad[d];
    }
    
    return rbf_->dd_basis(shape * dist,
                          shape * shape * grad2 * grad2,
                          shape * lap);
}

void RBF_Function::
output(pugi::xml_node &output_node) const
{
    pugi::xml_node rbf_node = output_node.append_child("rbf_function");
    
    XML_Functions::append_child(rbf_node, "standard", "function_type");
    XML_Functions::append_child(rbf_node, rbf_->description(), "rbf_type");
    XML_Functions::append_child(rbf_node, distance_->description(), "distance_type");
}

void RBF_Function::
check_class_invariants() const
{
    Assert(rbf_);
    Assert(distance_);
}

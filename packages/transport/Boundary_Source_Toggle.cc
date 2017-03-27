#include "Boundary_Source_Toggle.hh"

#include "XML_Node.hh"

using namespace std;

Boundary_Source_Toggle::
Boundary_Source_Toggle(bool local_include_boundary_source,
                       shared_ptr<Sweep_Operator> sweep):
    Sweep_Operator(sweep->sweep_type(),
                   sweep->transport_discretization()),
    local_include_boundary_source_(local_include_boundary_source),
    sweep_(sweep)
{
}

void Boundary_Source_Toggle::
apply(vector<double> &x) const
{
    sweep_->set_include_boundary_source(local_include_boundary_source_);
    (*sweep_)(x);
}

void Boundary_Source_Toggle::
output(XML_Node output_node) const
{
    sweep_->output(output_node);
}

#include "Source_Iteration.hh"

#include "Angular_Discretization.hh"
#include "Convergence_Measure.hh"
#include "Energy_Discretization.hh"
#include "Spatial_Discretization.hh"
#include "Transport_Discretization.hh"
#include "Vector_Operator.hh"
#include "XML_Node.hh"

using namespace std;

Source_Iteration::
Source_Iteration(Options options,
                 shared_ptr<Spatial_Discretization> spatial_discretization,
                 shared_ptr<Angular_Discretization> angular_discretization,
                 shared_ptr<Energy_Discretization> energy_discretization,
                 shared_ptr<Transport_Discretization> transport_discretization,
                 shared_ptr<Convergence_Measure> convergence,
                 shared_ptr<Vector_Operator> source_operator,
                 shared_ptr<Vector_Operator> flux_operator,
                 vector<shared_ptr<Vector_Operator> > value_operators):
    Solver(options.solver_print,
           Solver::Type::STEADY_STATE),
    options_(options),
    spatial_discretization_(spatial_discretization),
    angular_discretization_(angular_discretization),
    energy_discretization_(energy_discretization),
    transport_discretization_(transport_discretization),
    convergence_(convergence),
    source_operator_(source_operator),
    flux_operator_(flux_operator),
    value_operators_(value_operators)
{
    convergence_->set_tolerance(options.tolerance);
}

void Source_Iteration::
solve()
{
    int phi_size = transport_discretization_->phi_size();
    int number_of_augments = transport_discretization_->number_of_augments();

    // Initialize result
    result_ = make_shared<Result>();
    
    // Calculate first-flight source
    vector<double> q(phi_size + number_of_augments, 0);
    print_name("Initial source iteration");
    if (transport_discretization_->has_reflection())
    {
        double error = 1;
        double error_old = 1;
        vector<double> q_old;
        for (int it = 0; it < options_.max_source_iterations; ++it)
        {
            print_iteration(it);
            
            // Perform sweep to get new phi
            q_old = q;
            (*source_operator_)(q);
            
            // Get error
            error_old = error;
            error = convergence_->error(q,
                                        q_old);
            print_error(error);

            // Check convergence
            bool converged = convergence_->check(error,
                                                 error_old);
            if (converged)
            {
                result_->source_iterations = it + 1;
                print_convergence();
                break;
            }
        }
    }
    else
    {
        // Without reflection, only one application of operator is needed
        print_iteration(0);
        (*source_operator_)(q); 
        print_error(0);
        print_convergence();
   }

    // Zero out augments of first-flight source
    for (int i = phi_size; i < phi_size + number_of_augments; ++i)
    {
        q[i] = 0;
    }

    // Perform source iterations
    print_name("Source iteration");
    vector<double> x(q);
    double error = 1;
    {
        vector<double> x_old;
        double error_old = 1;
        for (int it = 0; it < options_.max_iterations; ++it)
        {
            print_iteration(it);

            // Perform sweep to get new phi
            x_old = x;
            (*flux_operator_)(x);

            // Add first-flight source to phi
            for (int i = 0; i < phi_size; ++i)
            {
                x[i] += q[i];
            }
            
            // Get error
            error_old = error;
            error = convergence_->error(x,
                                        x_old);;
            print_error(error);
            
            // Check convergence
            bool converged = convergence_->check(error,
                                                 error_old);
            if (converged)
            {
                result_->total_iterations = it + 1;
                print_convergence();
                break;
            }
        }
    }
    // If total iterations has not been changed, the result did not converge
    if (result_->total_iterations == -1)
    {
        result_->total_iterations = options_.max_iterations;
        print_failure();
    }

    // Remove augments from result
    x.resize(phi_size);

    // Store coefficients
    result_->coefficients = x;

    // Get result
    int number_of_values = value_operators_.size();
    result_->phi.resize(number_of_values);
    for (int i = 0; i < number_of_values; ++i)
    {
        vector<double> &phi = result_->phi[i];
        phi = x;
        (*value_operators_[i])(phi);
    }
}

void Source_Iteration::
output(XML_Node output_node) const
{
    // Output options
    output_node.set_attribute(options_.max_source_iterations,
                              "max_source_iterations");
    output_node.set_attribute(options_.max_iterations,
                               "max_iterations");
    output_node.set_attribute(options_.solver_print,
                              "solver_print");
    output_node.set_attribute(options_.tolerance,
                              "tolerance");
    
    // Output results
    output_result(output_node,
                  result_);
}

void Source_Iteration::
check_class_invariants() const
{
    Assert(spatial_discretization_);
    Assert(angular_discretization_);
    Assert(energy_discretization_);
    Assert(transport_discretization_);
    Assert(convergence_);
    Assert(source_operator_);
    Assert(flux_operator_);
    for (std::shared_ptr<Vector_Operator> oper : value_operators_)
    {
        Assert(oper);
    }
}

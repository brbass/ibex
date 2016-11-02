#include "Source_Iteration.hh"

#include <cmath>

#include "Angular_Discretization.hh"
#include "Augmented_Operator.hh"
#include "Check.hh"
#include "Discrete_To_Moment.hh"
#include "Energy_Discretization.hh"
#include "Fission.hh"
#include "Identity_Operator.hh"
#include "Moment_To_Discrete.hh"
#include "Multiplicative_Operator.hh"
#include "Preconditioner.hh"
#include "Scattering.hh"
#include "Sweep_Operator.hh"
#include "Spatial_Discretization.hh"
#include "Transport_Discretization.hh"
#include "Vector_Functions.hh"
#include "Vector_Operator_Functions.hh"
#include "XML_Functions.hh"

using namespace std;

namespace vf = Vector_Functions;

Source_Iteration::
Source_Iteration(int max_iterations,
                 int solver_print,
                 double tolerance,
                 shared_ptr<Spatial_Discretization> spatial_discretization,
                 shared_ptr<Angular_Discretization> angular_discretization,
                 shared_ptr<Energy_Discretization> energy_discretization,
                 shared_ptr<Transport_Discretization> transport_discretization,
                 shared_ptr<Sweep_Operator> sweeper,
                 shared_ptr<Discrete_To_Moment> discrete_to_moment,
                 shared_ptr<Moment_To_Discrete> moment_to_discrete,
                 shared_ptr<Scattering> scattering,
                 shared_ptr<Fission> fission,
                 shared_ptr<Preconditioner> preconditioner):
    Solver(solver_print),
    max_iterations_(max_iterations),
    total_iterations_(0),
    source_iterations_(0),
    tolerance_(tolerance),
    spatial_discretization_(spatial_discretization),
    angular_discretization_(angular_discretization),
    energy_discretization_(energy_discretization),
    transport_discretization_(transport_discretization),
    sweeper_(sweeper),
    discrete_to_moment_(discrete_to_moment),
    moment_to_discrete_(moment_to_discrete),
    scattering_(scattering),
    fission_(fission),
    preconditioner_(preconditioner)
{
    if (preconditioner_ == shared_ptr<Preconditioner>())
    {
        preconditioned_ = false;
    }
    else
    {
        preconditioned_ = true;
    }
}

/*
  Apply phi(l+1) = D(Linv(M(S phi(l)))) + D(Linv(q))
*/
void Source_Iteration::
solve_steady_state(vector<double> &x)
{
    shared_ptr<Vector_Operator> SI = make_shared<Source_Iterator>(*this);
    shared_ptr<Vector_Operator> FI = make_shared<Flux_Iterator>(*this);

    int phi_size = transport_discretization_->phi_size();
    int number_of_augments = transport_discretization_->number_of_augments();
    
    vector<double> q(phi_size + number_of_augments, 0);
    
    if (transport_discretization_->has_reflection())
    {
        vector<double> q_old;
        
        print_name("Initial source iteration");
        
        for (int it = 0; it < max_iterations_; ++it)
        {
            print_iteration(it);
            
            q_old = q;
            
            (*SI)(q);

            double error;
            bool converged = check_phi_convergence(q, q_old, error);
            print_error(error);
            
            if (converged)
            {
                source_iterations_ = it + 1;
                
                print_convergence();
                
                break;
            }
        }
        for (int i = phi_size; i < phi_size + number_of_augments; ++i)
        {
            q[i] = 0;
        }
        if (source_iterations_ == max_iterations_)
        {
            print_failure();
        }
    }
    else
    {
        (*SI)(q);
        
        source_iterations_ = 1;
    }
    
    x.resize(phi_size + number_of_augments, 0);
    vector<double> x_old;
    
    print_name("Source iteration");
    
    for (int it = 0; it < max_iterations_; ++it)
    {
        print_iteration(it);

        x_old = x;
        
        (*FI)(x);
        
        for (int i = 0; i < phi_size; ++i)
        {
            x[i] += q[i]; // add source
        }

        double error;
        bool converged = check_phi_convergence(x, x_old, error);
        print_error(error);
        
        if (converged)
        {
            total_iterations_ = it + 1;

            print_convergence();
            
            break;
        }
    }
    if (total_iterations_ == max_iterations_)
    {
        print_failure();
    }
    
    x.resize(phi_size); // remove augments
}

/*
  Check convergence of pointwise relative error in scalar flux
*/
bool Source_Iteration::
check_phi_convergence(vector<double> const &x, 
                      vector<double> const &x_old,
                      double &error) const
{
    int number_of_cells = spatial_discretization_->number_of_points();
    int number_of_nodes = spatial_discretization_->number_of_nodes();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_moments = angular_discretization_->number_of_moments();
    
    {
        int m = 0;
        
        for (int i = 0; i < number_of_cells; ++i)
        {
            for (int g = 0; g < number_of_groups; ++g)
            {
                for (int n = 0; n < number_of_nodes; ++n)
                {
                    int k = n + number_of_nodes * (g + number_of_groups * (m + number_of_moments * i));
                    error = abs(x[k] - x_old[k]) / (abs(x_old[k]) + tolerance_ * tolerance_);
                    
                    if (error > tolerance_)
                    {
                        return false;
                    }
                }
            }
        }
    }
    
    return true;
}

/*
  Check convergence of the k eigenvalue
*/
bool Source_Iteration::
check_k_convergence(double k,
                    double k_old,
                    double &error) const
{
    error = abs(k - k_old) / (abs(k_old) + tolerance_ * tolerance_);
    
    if (error > tolerance_)
    {
        return false;
    }
    else
    {
        return true;
    }
}

void Source_Iteration::
solve_k_eigenvalue(double &k_eigenvalue, 
                   vector<double> &x)
{
    int phi_size = transport_discretization_->phi_size();
    int number_of_augments = transport_discretization_->number_of_augments();

    k_eigenvalue = 1.0;
    x.resize(phi_size + number_of_augments, 1.0);

    shared_ptr<Vector_Operator> D
        = make_shared<Augmented_Operator>(number_of_augments,
                                          discrete_to_moment_);
    shared_ptr<Vector_Operator> M
        = make_shared<Augmented_Operator>(number_of_augments,
                                          moment_to_discrete_);
    shared_ptr<Vector_Operator> S
        = make_shared<Augmented_Operator>(number_of_augments,
                                          scattering_);
    shared_ptr<Vector_Operator> F
        = make_shared<Augmented_Operator>(number_of_augments,
                                          fission_,
                                          true); // zero out augments
    shared_ptr<Multiplicative_Operator> U
        = make_shared<Multiplicative_Operator>(phi_size + number_of_augments,
                                               k_eigenvalue);
    shared_ptr<Sweep_Operator> Linv = sweeper_;
    Assert(Linv);
    Linv->set_include_boundary_source(false);
    
    shared_ptr<Vector_Operator> T;
    switch (Linv->sweep_type())
    {
    case Sweep_Operator::Sweep_Type::MOMENT:
        T = Linv;
        break;
    case Sweep_Operator::Sweep_Type::ORDINATE:
        T = D * Linv * M;
        break;
    }
    
    shared_ptr<Vector_Operator> Op = T * (S + U * F);
    
    double k_eigenvalue_old;
    vector<double> x_old;
    
    print_name("Fixed point iteration");
    
    for (int it = 0; it < max_iterations_; ++it)
    {
        print_iteration(it);
        
        k_eigenvalue_old = k_eigenvalue;
        x_old = x;
        
        U->set_scalar(1. / k_eigenvalue);
        
        (*Op)(x);
        
        k_eigenvalue = update_k(k_eigenvalue_old,
                                x,
                                x_old);
        
        double error;
        bool converged = check_k_convergence(k_eigenvalue, k_eigenvalue_old, error);
        print_value(k_eigenvalue);
        print_error(error);
        
        if (converged)
        {
            total_iterations_ = it + 1;

            print_convergence();
            print_eigenvalue(k_eigenvalue);
            
            break;
        }
    }
    if (total_iterations_ == max_iterations_)
    {
        print_failure();
    }
    
    x.resize(phi_size); // remove augments
}

void Source_Iteration::
solve_time_dependent(vector<double> &x)
{
    AssertMsg(false, "not implemented");
}

double Source_Iteration::
update_k(double k_eigenvalue_old,
         vector<double> const &x,
         vector<double> const &x_old) const
{
    int phi_size = transport_discretization_->phi_size();

    vector<double> fis(x.begin(), x.begin() + phi_size);
    vector<double> fis_old(x_old.begin(), x_old.begin() + phi_size);
    vector<double> sca(x.begin(), x.begin() + phi_size);
    vector<double> sca_old(x_old.begin(), x_old.begin() + phi_size);
    
    (*fission_)(fis);
    (*fission_)(fis_old);
    (*scattering_)(sca);
    (*scattering_)(sca_old);

    double num = vf::magnitude(fis);
    // double den = vf::magnitude(fis_old) / k_eigenvalue_old;
    double den = vf::magnitude(vf::subtract(vf::multiply(fis_old,
                                                         1. / k_eigenvalue_old),
                                            vf::subtract(sca,
                                                         sca_old)));
    
    return num / den;
}

void Source_Iteration::
output(pugi::xml_node &output_node) const
{
    pugi::xml_node source = output_node.append_child("solution_method");

    XML_Functions::append_child(source, "source_iteration", "type");
    XML_Functions::append_child(source, max_iterations_, "max_iterations");
    XML_Functions::append_child(source, source_iterations_, "source_iterations");
    XML_Functions::append_child(source, total_iterations_, "total_iterations");
    XML_Functions::append_child(source, tolerance_, "tolerance");
}

Source_Iteration::Source_Iterator::
Source_Iterator(Source_Iteration const &si):
    Vector_Operator(si.transport_discretization_->phi_size() + si.transport_discretization_->number_of_augments(),
                    si.transport_discretization_->phi_size() + si.transport_discretization_->number_of_augments()),
    si_(si)
{
}

void Source_Iteration::Source_Iterator::
apply(vector<double> &x) const
{
    int phi_size = si_.transport_discretization_->phi_size();
    int number_of_augments = si_.transport_discretization_->number_of_augments();
    
    shared_ptr<Vector_Operator> D = make_shared<Augmented_Operator>(number_of_augments, si_.discrete_to_moment_);
    shared_ptr<Vector_Operator> M = make_shared<Augmented_Operator>(number_of_augments, si_.moment_to_discrete_);
    
    shared_ptr<Sweep_Operator> Linv = si_.sweeper_;
    Linv->set_include_boundary_source(true);
    
    shared_ptr<Vector_Operator> T;
    switch (Linv->sweep_type())
    {
    case Sweep_Operator::Sweep_Type::MOMENT:
        T = Linv;
        break;
    case Sweep_Operator::Sweep_Type::ORDINATE:
        T = D * Linv * M;
        break;
    }

    shared_ptr<Vector_Operator> E;
    if (si_.preconditioned_)
    {
        shared_ptr<Preconditioner> Cinv =si_.preconditioner_;
        shared_ptr<Sweep_Operator> sweeper = Cinv->sweeper();
        sweeper->set_include_boundary_source(false);
        switch (sweeper->sweep_type())
        {
        case Sweep_Operator::Sweep_Type::MOMENT:
            E = Cinv;
            break;
        case Sweep_Operator::Sweep_Type::ORDINATE:
            E = D * Cinv * M;
            break;
        }
    }

    int number_of_points = si_.spatial_discretization_->number_of_points();
    int number_of_nodes = si_.spatial_discretization_->number_of_nodes();
    int number_of_groups = si_.energy_discretization_->number_of_groups();
    int number_of_moments = si_.angular_discretization_->number_of_moments();
    
    for (int i = 0; i < number_of_points; ++i)
    {
        shared_ptr<Material> material = si_.spatial_discretization_->point(i)->material();
        vector<double> const internal_source = material->internal_source();
        
        for (int g = 0; g < number_of_groups; ++g)
        {
            for (int n = 0; n < number_of_nodes; ++n)
            {
                {
                    int m = 0;

                    int k = n + number_of_nodes * (g + number_of_groups * (m + number_of_moments * i));

                    x[k] = internal_source[g];
                }
                
                for (int m = 1; m < number_of_moments; ++m)
                {
                    int k = n + number_of_nodes * (g + number_of_groups * (m + number_of_moments * i));
                    
                    x[k] = 0;
                }
            }
        }
    }
    
    shared_ptr<Vector_Operator> Op;

    if (si_.preconditioned_)
    {
        shared_ptr<Vector_Operator> S
            = make_shared<Augmented_Operator>(number_of_augments,
                                              si_.scattering_);
        shared_ptr<Vector_Operator> F
            = make_shared<Augmented_Operator>(number_of_augments,
                                              si_.fission_,
                                              true); // Only one of S or F needs to include the augments
        shared_ptr<Vector_Operator> I
            = make_shared<Identity_Operator>(phi_size + number_of_augments);
        
        Op = (I + E * (S + F)) * T;
    }
    else
    {
        Op = T;
    }
    
    (*Op)(x);
}

Source_Iteration::Flux_Iterator::
Flux_Iterator(Source_Iteration const &si):
    Vector_Operator(si.transport_discretization_->phi_size() + si.transport_discretization_->number_of_augments(),
                    si.transport_discretization_->phi_size() + si.transport_discretization_->number_of_augments()),
    si_(si)
{
    
}

void Source_Iteration::Flux_Iterator::
apply(vector<double> &x) const
{
    int phi_size = si_.transport_discretization_->phi_size();
    int number_of_augments = si_.transport_discretization_->number_of_augments();

    shared_ptr<Vector_Operator> D
        = make_shared<Augmented_Operator>(number_of_augments,
                                          si_.discrete_to_moment_);
    shared_ptr<Vector_Operator> M
        = make_shared<Augmented_Operator>(number_of_augments,
                                          si_.moment_to_discrete_);

    shared_ptr<Sweep_Operator> Linv = si_.sweeper_;
    Linv->set_include_boundary_source(false);

    shared_ptr<Vector_Operator> T;
    switch (Linv->sweep_type())
    {
    case Sweep_Operator::Sweep_Type::MOMENT:
        T = Linv;
        break;
    case Sweep_Operator::Sweep_Type::ORDINATE:
        T = D * Linv * M;
        break;
    }

    shared_ptr<Vector_Operator> E;
    if (si_.preconditioned_)
    {
        shared_ptr<Preconditioner> Cinv = si_.preconditioner_;
        shared_ptr<Sweep_Operator> sweeper = Cinv->sweeper();
        sweeper->set_include_boundary_source(false);
        switch (sweeper->sweep_type())
        {
        case Sweep_Operator::Sweep_Type::MOMENT:
            E = Cinv;
            break;
        case Sweep_Operator::Sweep_Type::ORDINATE:
            E = D * Cinv * M;
            break;
        }
    }
    
    shared_ptr<Vector_Operator> S
        = make_shared<Augmented_Operator>(number_of_augments,
                                          si_.scattering_);
    shared_ptr<Vector_Operator> F
        = make_shared<Augmented_Operator>(number_of_augments,
                                          si_.fission_,
                                          true); // Only one of S or F needs to include the augments
    shared_ptr<Vector_Operator> I
        = make_shared<Identity_Operator>(phi_size + number_of_augments);

    shared_ptr<Vector_Operator> Op;
    
    if (si_.preconditioned_)
    {
        Op = I + (I + E * (S + F)) * (T * (S + F) - I);
    }
    else
    {
        Op = T * (S + F);
    }
    
    (*Op)(x);
}


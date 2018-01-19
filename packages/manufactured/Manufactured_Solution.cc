#include "Manufactured_Solution.hh"

#include "Angular_Discretization.hh"
#include "Angular_Discretization_Factory.hh"
#include "Check.hh"
#include "Energy_Discretization.hh"

using namespace std;

Manufactured_Solution::
Manufactured_Solution(shared_ptr<Angular_Discretization> angular,
                      shared_ptr<Energy_Discretization> energy):
    angular_(angular),
    energy_(energy)
{
    Assert(angular);
    Assert(energy);
    
    // Get temporary angular discretization for streaming coefficients
    int dimension = angular_->dimension();
    int number_of_scattering_moments = angular_->number_of_scattering_moments();
    Angular_Discretization_Factory factory;
    int angular_rule = (dimension == 1 ? 256 : 6);
    shared_ptr<Angular_Discretization> coeff_angular
        = factory.get_angular_discretization(dimension,
                                             number_of_scattering_moments,
                                             angular_rule);

    // Get streaming coefficients
    coeff_angular->manufactured_coefficients(streaming_size_,
                                             streaming_indices_,
                                             streaming_coefficients_);
}

vector<double> Manufactured_Solution::
get_source(vector<double> const &solution,
           vector<double> const &grad_solution,
           vector<double> const &sigma_t,
           vector<double> const &sigma_s) const
{
    // Get size and indexing information
    int dimension = angular_->dimension();
    int number_of_moments = angular_->number_of_moments();
    int number_of_scattering_moments = angular_->number_of_scattering_moments();
    int number_of_groups = energy_->number_of_groups();
    vector<int> const scattering_indices = angular_->scattering_indices();

    // Check sizes
    Check(solution.size() == number_of_moments * number_of_groups);
    Check(grad_solution.size() == number_of_moments * number_of_groups * dimension);
    Check(sigma_t.size() == number_of_groups);
    Check(sigma_s.size() == number_of_groups * number_of_groups * number_of_scattering_moments);

    // Initialize source to zero
    vector<double> source(number_of_groups * number_of_moments, 0);

    // Calculate source from group and moment gf and mf to group and moment gt and mt
    for (int mt = 0; mt < number_of_moments; ++mt)
    {
        // Get scattering index
        int l = scattering_indices[mt];
        
        // Get coefficient information for this moment
        int const &local_size = streaming_size_[mt];
        vector<int> const &local_indices = streaming_indices_[mt];
        vector<double> const &local_coefficients = streaming_coefficients_[mt];
        
        for (int gt = 0; gt < number_of_groups; ++ gt)
        {
            // Get element of source
            int k_src = gt + number_of_groups * mt;
            double &sum = source[k_src];

            // Add streaming term
            for (int s = 0; s < local_size; ++s)
            {
                int mf = local_indices[s];
                
                for (int d = 0; d < dimension; ++d)
                {
                    int k_coeff = d + dimension * s;
                    int k_sol = d + dimension * (gt + number_of_groups * mf);
                    sum += local_coefficients[k_coeff] * grad_solution[k_sol];
                }
            }
            
            // Add collision term
            {
                int k_sol = gt + number_of_groups * mt;
                sum += sigma_t[gt] * solution[k_sol];
            }
            
            // Subtract combined scattering and fission source
            for (int gf = 0; gf < number_of_groups; ++gf)
            {
                int k_sol = gf + number_of_groups * mt;
                int k_ss = gf + number_of_groups * (gt + number_of_groups * l);
                sum -= sigma_s[k_ss] * solution[k_sol];
            }
        }
    }
    
    return source;
}


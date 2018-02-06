#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <mpi.h>

#include "Heat_Transfer_Data.hh"
#include "Heat_Transfer_Factory.hh"
#include "Heat_Transfer_Integration.hh"
#include "Heat_Transfer_Solution.hh"
#include "Heat_Transfer_Solve.hh"
#include "Weak_Spatial_Discretization.hh"

using namespace std;

// Local class for given heat transfer data
class Constant_Heat_Transfer_Data : public Heat_Transfer_Data
{
public:
    
    Constant_Heat_Transfer_Data(double r1,
                                double r2,
                                double k1,
                                double k2,
                                double h,
                                double q1,
                                double q2,
                                double tinf):
        Heat_Transfer_Data(),
        r1_(r1),
        r2_(r2),
        k1_(k1),
        k2_(k2),
        h_(h),
        q1_(q1),
        q2_(q2),
        tinf_(tinf)
    {
    }
    virtual double conduction(std::vector<double> const &position) const override
    {
        return position[0] < r1_ ? k1_ : k2_;
    }
    virtual double convection(std::vector<double> const &position) const override
    {
        return h_;
    }
    virtual double source(std::vector<double> const &position) const override
    {
        return position[0] < r1_ ? q1_ : q2_;
    }
    virtual double temperature_inf(std::vector<double> const &position) const override
    {
        return tinf_;
    }

    double solution(std::vector<double> const &position) const
    {
        double r = position[0];
        if (r < r1_)
        {
            return (2*k1_*k2_*q1_*pow(r1_,2) - 2*k1_*k2_*q2_*pow(r1_,2) - h_*k2_*q1_*pow(r,2)*r2_ + h_*k2_*q1_*pow(r1_,2)*r2_ - h_*k1_*q2_*pow(r1_,2)*r2_ + 2*k1_*k2_*q2_*pow(r2_,2) + h_*k1_*q2_*pow(r2_,3) + 4*h_*k1_*k2_*r2_*tinf_ + 2*h_*k1_*(-q1_ + q2_)*pow(r1_,2)*r2_*log(r1_) + 2*h_*k1_*(q1_ - q2_)*pow(r1_,2)*r2_*log(r2_))/(4.*h_*k1_*k2_*r2_);
        }
        else
        {
            return (2*k2_*q1_*pow(r1_,2) - 2*k2_*q2_*pow(r1_,2) - h_*q2_*pow(r,2)*r2_ + 2*k2_*q2_*pow(r2_,2) + h_*q2_*pow(r2_,3) + 4*h_*k2_*r2_*tinf_ + 2*h_*(-q1_ + q2_)*pow(r1_,2)*r2_*log(r) + 2*h_*(q1_ - q2_)*pow(r1_,2)*r2_*log(r2_))/(4.*h_*k2_*r2_);
        }
    }
    
private:
    
    double r1_;
    double r2_;
    double k1_;
    double k2_;
    double h_;
    double q1_;
    double q2_;
    double tinf_;
};

int test_constant(int number_of_points,
                  double radius_num_intervals,
                  double length1,
                  double length2,
                  double conduction1,
                  double conduction2,
                  double convection,
                  double source1,
                  double source2,
                  double temperature_inf)
{
    int checksum = 0;
    
    // Get data
    shared_ptr<Constant_Heat_Transfer_Data> data
        = make_shared<Constant_Heat_Transfer_Data>(length1,
                                                   length2,
                                                   conduction1,
                                                   conduction2,
                                                   convection,
                                                   source1,
                                                   source2,
                                                   temperature_inf);
    
    // Get solver
    Heat_Transfer_Factory factory;
    shared_ptr<Heat_Transfer_Solve> solver
        = factory.get_cylindrical_1d(number_of_points,
                                     radius_num_intervals,
                                     length2,
                                     true, // basis mls
                                     true, // weight mls
                                     "wendland12",
                                     "wendland12",
                                     data);
    
    // Solve problem
    shared_ptr<Heat_Transfer_Solution> solution
        = solver->solve();
    
    vector<double> test_position(1, 0.5);
    cout << solution->solution(test_position) << endl;
    cout << data->solution(test_position) << endl;

    return checksum;
}
                  

int main(int argc, char **argv)
{
    int checksum = 0;
    
    MPI_Init(&argc, &argv);

    int number_of_points = 100;
    double radius_num_intervals = 3;
    double length1 = 1;
    double length2 = 1;
    double conduction1 = 0.5;
    double conduction2 = 1.0;
    double convection = 2.0;
    double source1 = 1.0;
    double source2 = 2.0;
    double temperature_inf = 500;
    checksum += test_constant(number_of_points,
                              radius_num_intervals,
                              length1,
                              length2,
                              conduction1,
                              conduction2,
                              convection,
                              source1,
                              source2,
                              temperature_inf);
    
    MPI_Finalize();
    
    return checksum;
}


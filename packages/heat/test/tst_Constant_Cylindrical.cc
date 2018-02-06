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
    
    Constant_Heat_Transfer_Data(double length1,
                                double length2,
                                double conduction1,
                                double conduction2,
                                double convection,
                                double source1,
                                double source2,
                                double temperature_inf):
        Heat_Transfer_Data(),
        length1_(length1),
        length2_(length2),
        conduction1_(conduction1),
        conduction2_(conduction2),
        convection_(convection),
        source1_(source1),
        source2_(source2),
        temperature_inf_(temperature_inf)
    {
    }
    virtual double conduction(std::vector<double> const &position) const override
    {
        return position[0] < length1_ ? conduction1_ : conduction2_;
    }
    virtual double convection(std::vector<double> const &position) const override
    {
        return convection_;
    }
    virtual double source(std::vector<double> const &position) const override
    {
        return position[0] < length1_ ? source1_ : source2_;
    }
    virtual double temperature_inf(std::vector<double> const &position) const override
    {
        return temperature_inf_;
    }

private:
    
    double length1_;
    double length2_;
    double conduction1_;
    double conduction2_;
    double convection_;
    double source1_;
    double source2_;
    double temperature_inf_;
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
    shared_ptr<Heat_Transfer_Data> data
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


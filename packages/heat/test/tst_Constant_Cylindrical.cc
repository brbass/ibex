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
#include "Weak_Spatial_Discretization_Parser.hh"
#include "XML_Document.hh"
#include "XML_Node.hh"

using namespace std;

// Local class for given heat transfer data
class Constant_Heat_Transfer_Data : public Heat_Transfer_Data
{
public:
    
    Constant_Heat_Transfer_Data(int dimension,
                                double r1,
                                double r2,
                                double k1,
                                double k2,
                                double h,
                                double q1,
                                double q2,
                                double tinf):
        Heat_Transfer_Data(),
        dimension_(dimension),
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
        return radius2(position) < r1_ * r1_ ? k1_ : k2_;
    }
    virtual double convection(std::vector<double> const &position) const override
    {
        return h_;
    }
    virtual double source(std::vector<double> const &position) const override
    {
        return radius2(position) < r1_ * r1_ ? q1_ : q2_;
    }
    virtual double temperature_inf(std::vector<double> const &position) const override
    {
        return tinf_;
    }

    double solution(std::vector<double> const &position) const
    {
        double r = sqrt(radius2(position));
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

    double radius2(std::vector<double> const &position) const
    {
        switch (dimension_)
        {
        case 1:
            return position[0] * position[0];
        case 2:
            return position[0] * position[0] + position[1] * position[1];
        }
    }
    
    int dimension_;
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
        = make_shared<Constant_Heat_Transfer_Data>(1, // dimension
                                                   length1,
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
                                     "wendland11",
                                     "wendland11",
                                     data);
    
    // Solve problem
    shared_ptr<Heat_Transfer_Solution> solution
        = solver->solve();
    
    vector<double> test_position(1, 0.046);
    cout << solution->solution(test_position) << endl;
    cout << data->solution(test_position) << endl;

    return checksum;
}

int test_constant_2d(XML_Node input_node,
                     double length1,
                     double length2,
                     double conduction1,
                     double conduction2,
                     double convection,
                     double source1,
                     double source2,
                     double temperature_inf)
{
    // Get solid geometry
    int dimension = 2;
    vector<vector<double> > limits = {{-length2, length2},
                                      {-length2, length2}};
    shared_ptr<Solid_Geometry> solid;
    vector<shared_ptr<Cartesian_Plane> > surfaces;
    Heat_Transfer_Factory factory;
    factory.get_solid(dimension,
                      limits,
                      solid,
                      surfaces);
    
    // Get weak spatial discretization
    Weak_Spatial_Discretization_Parser spatial_parser(solid,
                                                      surfaces);
    shared_ptr<Weak_Spatial_Discretization> spatial
        = spatial_parser.get_weak_discretization(input_node.get_child("spatial_discretization"));
    
    // Get data
    shared_ptr<Constant_Heat_Transfer_Data> data
        = make_shared<Constant_Heat_Transfer_Data>(dimension,
                                                   length1,
                                                   length2,
                                                   conduction1,
                                                   conduction2,
                                                   convection,
                                                   source1,
                                                   source2,
                                                   temperature_inf);
    
    // Get heat transfer integration options
    shared_ptr<Heat_Transfer_Integration_Options> integration_options
        = make_shared<Heat_Transfer_Integration_Options>();
    integration_options->geometry = Heat_Transfer_Integration_Options::Geometry::CYLINDRICAL_2D;
    
    // Perform integration
    shared_ptr<Heat_Transfer_Integration> integration
        = make_shared<Heat_Transfer_Integration>(integration_options,
                                                 data,
                                                 spatial);
    
    // Initialize heat transfer solver
    shared_ptr<Heat_Transfer_Solve> solver
        = make_shared<Heat_Transfer_Solve>(integration,
                                           spatial);
    
    // Solve problem
    shared_ptr<Heat_Transfer_Solution> solution
        = solver->solve();

    // Check a spot
    vector<double> test_position(dimension, 0.046);
    test_position[1] = 0;
    cout << solution->solution(test_position) << endl;
    cout << data->solution(test_position) << endl;
}

int main(int argc, char **argv)
{
    int checksum = 0;
    
    MPI_Init(&argc, &argv);

    if (argc == 1)
    {
        int number_of_points = 400;
        double radius_num_intervals = 3.1;
        double length1 = 1.4;
        double length2 = 2;
        double conduction1 = 0.0007;
        double conduction2 = 0.05;
        double convection = 3.;
        double source1 = 3.;
        double source2 = 0.;
        double temperature_inf = 600;
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
    }
    else
    {
        Assert(argc == 2);
        string input_filename = argv[1];
        XML_Document input_file(input_filename);
        XML_Node input_node = input_file.get_child("input");
        
        double length1 = 1.4;
        double length2 = 2;
        double conduction1 = 0.0007;
        double conduction2 = 0.05;
        double convection = 3.;
        double source1 = 3.;
        double source2 = 0.;
        double temperature_inf = 600;
        checksum += test_constant_2d(input_node,
                                     length1,
                                     length2,
                                     conduction1,
                                     conduction2,
                                     convection,
                                     source1,
                                     source2,
                                     temperature_inf);
    }
    MPI_Finalize();
    
    return checksum;
}


#include <iostream>

#include "Cartesian_Integration_1D.hh"
#include "Cartesian_Integration_2D.hh"
#include "Cartesian_Integration_3D.hh"
#include "Check_Equality.hh"
#include "Cylindrical_Integration_2D.hh"
#include "Random_Number_Generator.hh"
#include "Spherical_Integration_3D.hh"

namespace ce = Check_Equality;

using namespace std;

namespace
{
    Random_Number_Generator<double> rng(-1,   // min
                                        1,    // max
                                        482); // seed
}

int test_cartesian_1d(int order, int num_tests)
{
    int checksum = 0;

    double tolerance = 1e-14;

    Cartesian_Integration_1D integration(order);

    // Check volume
    {
        Integrand_1D func = [=](double x){return 1;};
        
        for (int i = 0; i < num_tests; ++i)
        {
            double x1 = rng.scalar();
            double x2 = rng.scalar() + x1;
            
            double integral = integration.integrate(func,
                                                    x1,
                                                    x2);
            double analytic = x2 - x1;
        
            if (!ce::approx(integral, analytic, tolerance))
            {
                cerr << "cart integration 1d volume incorrect" << endl;
                cerr << "\texpected: " << analytic << "\tcalculated: " << integral << endl;
                checksum += 1;
            }
        }
    }
    
    // Check quartic functions
    {
        for (int i = 0; i < num_tests; ++i)
        {
            vector<double> a = rng.vector(5);
            
            Integrand_1D func = [=](double x){return a[0] + a[1]*x + a[2]*x*x + a[3]*x*x*x + a[4]*x*x*x*x;};
            
            double x1 = rng.scalar();
            double x2 = rng.scalar() + x1;
            
            double integral = integration.integrate(func,
                                                    x1,
                                                    x2);
            double analytic = 1./60.*(-60*x1*a[0]-30*x1*x1*a[1]-20*x1*x1*x1*a[2]-15*x1*x1*x1*x1*a[3]+5*x2*(12*a[0]+x2*(6*a[1]+4*x2*a[2]+3*x2*x2*a[3]))-12*x1*x1*x1*x1*x1*a[4]+12*x2*x2*x2*x2*x2*a[4]);
            
            if (!ce::approx(integral, analytic, tolerance))
            {
                cerr << "cart integration 1d integral incorrect, index: " << i << endl;
                cerr << "\texpected: " << analytic << "\tcalculated: " << integral << endl;
                checksum += 1;
            }
        }
    }
    
    return checksum;
}

int test_cartesian_2d(int order, int num_tests)
{
    int checksum = 0;
    
    double tolerance = 1e-12;
    
    Cartesian_Integration_2D integration(order);

    // Check volume
    {
        Integrand_2D func = [=](double x, double y){return 1;};

        for (int i = 0; i < num_tests; ++i)
        {
            double x1 = rng.scalar();
            double x2 = rng.scalar() + x1;
            double y1 = rng.scalar();
            double y2 = rng.scalar() + y1;
            cout << x2 << "\t" << y2 << endl;
            double integral = integration.integrate(func,
                                                    x1,
                                                    x2,
                                                    y1,
                                                    y2);
            double analytic = (x2 - x1) * (y2 - y1);
            
            if (!ce::approx(integral, analytic, tolerance))
            {
                cerr << "cart integration 2d volume incorrect, index: " << i << endl;
                cerr << "\texpected: " << analytic << "\tcalculated: " << integral << endl;
                checksum += 1;
            }
        }
    }
    
    // Check quadratic functions
    {
        for (int i = 0; i < num_tests; ++i)
        {
            vector<double> a = rng.vector(5);
            
            Integrand_2D func = [=](double x, double y){return a[0] + a[1] * x + a[2] * y + a[3] * x * x + a[4] * y * y + a[5] * x * y;};
            
            double x1 = rng.scalar();
            double x2 = rng.scalar() + x1;
            double y1 = rng.scalar();
            double y2 = rng.scalar() + y1;
            
            double integral = integration.integrate(func,
                                                    x1,
                                                    x2,
                                                    y1,
                                                    y2);
            double analytic = 1./12.*(x1-x2)*(y1-y2)*(12*a[0]+6*(x1+x2)*a[1]+6*(y1+y2)*a[2]+4*(x1*x1+x1*x2+x2*x2)*a[3]+4*(y1*y1+y1*y2+y2*y2)*a[4]+3*(x1+x2)*(y1+y2)*a[5]);
            
            if (!ce::approx(integral, analytic, tolerance))
            {
                cerr << "cart integration 2d integral incorrect, index: " << i << endl;
                cerr << "\texpected: " << analytic << "\tcalculated: " << integral << endl;
                checksum += 1;
            }
        }
    }
    
    return checksum;
}

int test_cylindrical_2d(int order, int num_tests)
{
    int checksum = 0;
    
    double tolerance = 1e-12;

    Cylindrical_Integration_2D integration(order);

    // Check volume
    {
        Integrand_2D func = [=](double x, double y){return 1;};

        for (int i = 0; i < num_tests; ++i)
        {
            
            double x0 = rng.scalar();
            double y0 = rng.scalar();
            double r1 = rng.scalar();
        
            double integral = integration.integrate(func,
                                                    x0,
                                                    y0,
                                                    r1);
            double analytic = M_PI * r1 * r1;
            
            if (!ce::approx(integral, analytic, tolerance))
            {
                cerr << "cyl integration 2d volume incorrect, index: " << i << endl;
                cerr << "\texpected: " << analytic << "\tcalculated: " << integral << endl;
                checksum += 1;
            }
        }
    }
    
    // Check Cartesian quadratic functions
    {
        for (int i = 0; i < num_tests; ++i)
        {
            vector<double> a = rng.vector(5);
            
            Integrand_2D func = [=](double x, double y){return a[0] + a[1] * x + a[2] * y + a[3] * x * x + a[4] * y * y + a[5] * x * y;};
            
            double x0 = rng.scalar();
            double y0 = rng.scalar();
            double r1 = rng.scalar();
            
            double integral = integration.integrate(func,
                                                    x0,
                                                    y0,
                                                    r1);
            double analytic = 0.25*M_PI*r1*r1*(4*a[0]+r1*r1*(a[3]+a[4])+4*y0*(a[2]+y0*a[4])+4*x0*(a[1]+x0*a[3]+y0*a[5]));
            
            if (!ce::approx(integral, analytic, tolerance))
            {
                cerr << "cylind integration 2d integral incorrect, index: " << i << endl;
                cerr << "\texpected: " << analytic << "\tcalculated: " << integral << endl;
                checksum += 1;
            }
        }
    }
    
    return checksum;
}

int test_cartesian_3d(int order, int num_tests)
{
    int checksum = 0;
    
    double tolerance = 1e-12;
    
    Cartesian_Integration_3D integration(order);
    
    // Check volume
    {
        Integrand_3D func = [=](double x, double y, double z){return 1;};
        
        for (int i = 0; i < num_tests; ++i)
        {
            double x1 = rng.scalar();
            double x2 = rng.scalar() + x1;
            double y1 = rng.scalar();
            double y2 = rng.scalar() + y1;
            double z1 = rng.scalar();
            double z2 = rng.scalar() + z1;
            
            double integral = integration.integrate(func,
                                                    x1,
                                                    x2,
                                                    y1,
                                                    y2,
                                                    z1,
                                                    z2);
            double analytic = (x2 - x1) * (y2 - y1) * (z2 - z1);
            
            if (!ce::approx(integral, analytic, tolerance))
            {
                cerr << "cart integration 3d volume incorrect, index: " << i << endl;
                cerr << "\texpected: " << analytic << "\tcalculated: " << integral << endl;
                checksum += 1;
            }
        }
    }
    
    // Check quadratic functions (accidentally left out the "z" polynomial)
    {
        for (int i = 0; i < num_tests; ++i)
        {
            vector<double> a = rng.vector(9);
            
            Integrand_3D func = [=](double x, double y, double z){return a[0] + a[1] * x + a[2] * y + a[3] * x * x + a[4] * y * y + a[5] * z * z + a[6] * x * y + a[7] * y * z + a[8] * x * z;};
            
            double x1 = rng.scalar();
            double x2 = rng.scalar() + x1;
            double y1 = rng.scalar();
            double y2 = rng.scalar() + y1;
            double z1 = rng.scalar();
            double z2 = rng.scalar() + z1;
            
            double integral = integration.integrate(func,
                                                    x1,
                                                    x2,
                                                    y1,
                                                    y2,
                                                    z1,
                                                    z2);
            double analytic = -1./12.*(x1-x2)*(y1-y2)*(z1-z2)*(12*a[0]+6*(y1+y2)*a[2]+4*x1*x1*a[3]+4*x2*x2*a[3]+4*(y1*y1+y1*y2+y2*y2)*a[4]+4*(z1*z1+z1*z2+z2*z2)*a[5]+3*(y1+y2)*(z1+z2)*a[7]+3*x2*(2*a[1]+(y1+y2)*a[6]+(z1+z2)*a[8])+x1*(6*a[1]+4*x2*a[3]+3*(y1+y2)*a[6]+3*(z1+z2)*a[8]));
            
            if (!ce::approx(integral, analytic, tolerance))
            {
                cerr << "cart integration 3d integral incorrect, index: " << i << endl;
                cerr << "\texpected: " << analytic << "\tcalculated: " << integral << endl;
                checksum += 1;
            }
        }
    }

    return checksum;
}

int test_spherical_3d(int order, int num_tests)
{
    int checksum = 0;
    
    double tolerance = 1e-12;
    
    Spherical_Integration_3D integration(order);
    
    // Check volume
    {
        Integrand_3D func = [=](double x, double y, double z){return 1;};
        
        for (int i = 0; i < num_tests; ++i)
        {
            double x0 = rng.scalar();
            double y0 = rng.scalar();
            double z0 = rng.scalar();
            double rmax = rng.scalar();
            
            double integral = integration.integrate(func,
                                                    x0,
                                                    y0,
                                                    z0,
                                                    rmax);
            double analytic = 4./3.*M_PI*rmax*rmax*rmax;
            
            if (!ce::approx(integral, analytic, tolerance))
            {
                cerr << "sph integration 3d volume incorrect, index: " << i << endl;
                cerr << "\texpected: " << analytic << "\tcalculated: " << integral << endl;
                checksum += 1;
            }
        }
    }
    
    // Check Cartesian quadratic functions
    {
        for (int i = 0; i < num_tests; ++i)
        {
            vector<double> a = rng.vector(9);
            
            Integrand_3D func = [=](double x, double y, double z){return a[0] + a[1] * x + a[2] * y + a[3] * x * x + a[4] * y * y + a[5] * z * z + a[6] * x * y + a[7] * y * z + a[8] * x * z;};
            
            double x0 = rng.scalar();
            double y0 = rng.scalar();
            double z0 = rng.scalar();
            double r1 = rng.scalar();
            
            double integral = integration.integrate(func,
                                                    x0,
                                                    y0,
                                                    z0,
                                                    r1);
            double analytic = 4./3.*M_PI*(1./5.*r1*r1*r1*r1*r1*(a[3]+a[4]+a[5])+r1*r1*r1*(a[0]+x0*x0*a[3]+z0*z0*a[5]+y0*(a[2]+y0*a[4]+z0*a[7])+x0*(a[1]+y0*a[6]+z0*a[8])));
            
            if (!ce::approx(integral, analytic, tolerance))
            {
                cerr << "sph integration cartesain 3d integral incorrect, index: " << i << endl;
                cerr << "\texpected: " << analytic << "\tcalculated: " << integral << endl;
                checksum += 1;
            }
        }
    }
    
    return checksum;
}

int main()
{
    int checksum = 0;
    
    int num_tests = 100;
    
    checksum += test_cartesian_1d(16, num_tests);
    // checksum += test_cartesian_2d(16, num_tests);
    checksum += test_cartesian_3d(16, num_tests);
    // checksum += test_cylindrical_2d(16, num_tests);
    checksum += test_spherical_3d(16, num_tests);
    
    return checksum;
}

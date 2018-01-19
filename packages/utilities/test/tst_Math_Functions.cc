#include <cmath>
#include <iostream>
#include <limits>

#include "Check_Equality.hh"
#include "Math_Functions.hh"
#include "Random_Number_Generator.hh"
#include "Vector_Functions.hh"

namespace ce = Check_Equality;
namespace mf = Math_Functions;
namespace vf = Vector_Functions;

using namespace std;

double tolerance = 10 * numeric_limits<double>::epsilon();

int test_factorial()
{
    int checksum = 0;

    int result = 479001600;
    if (!ce::equal(mf::factorial(12), result))
    {
        cout << "factorial failed";
        cout << endl;
        checksum += 1;
    }

    return checksum;
}

int test_legendre_p()
{
    int checksum = 0;

    double exp1 = 2565501./16.;
    double cal1 = mf::legendre_polynomial(12, sqrt(3));
    if (!ce::approx(cal1,
                    exp1,
                    exp1 * tolerance))
    {
        cout << "legendre_polynomial failed";
        cout << endl;
        cout << "\texpected: ";
        cout << exp1;
        cout << "\tactual: ";
        cout << cal1;
        cout << endl;
        checksum += 1;
    }

    double exp2 = 240625./243.;
    if (!ce::approx(mf::legendre_polynomial(7, 4, 2./3.),
                    exp2,
                    exp2 * tolerance))
    {
        cout << "associated legendre_polynomial failed";
        cout << endl;
        checksum += 1;
    }
    
    return checksum;
}

int test_spherical_h()
{
    int checksum = 0;

    // Test specific case
    {
        int const l = 3;
        int const m = -3;
        
        vector<double> const car = {4. / sqrt(23), 2. / sqrt(23), sqrt(3. / 23)};
        vector<double> const sph = {4 / sqrt(23.), acos(2. / sqrt(7.))};

        double const result = 9. * sqrt(15. / 46.) / 46.;
        
        // Check spherical and xyz transfer
        
        vector<double> new_car(3);
        mf::spherical_to_xyz(sph[0], sph[1],
                             new_car[0], new_car[1], new_car[2]);

        if (!ce::approx(car, new_car, tolerance))
        {
            cout << "spherical_to_xyz failed";
            cout << endl;
            checksum += 1;
        }
        
        vector<double> new_sph(2);
        mf::xyz_to_spherical(car[0], car[1], car[2],
                             new_sph[0], new_sph[1]);
        
        if (!ce::approx(sph, new_sph, tolerance))
        {
            cout << "xyz_to_spherical failed";
            cout << endl;
            checksum += 1;
        }
        
        // Check against known result
        
        double rec_sph = mf::spherical_harmonic_rec(l, m,
                                                    sph[0], sph[1]);


        if (!ce::approx(rec_sph, result, tolerance))
        {
            cout << "spherical recursive failed";
            cout << endl;
            checksum += 1;
        }
        
        double rec_car = mf::spherical_harmonic_rec(l, m,
                                                    car[0], car[1], car[2]);
        
        if (!ce::approx(rec_car, result, tolerance))
        {
            cout << "cartesian recursive failed";
            cout << endl;
            checksum += 1;
        }
        
        double tab_car = mf::spherical_harmonic(l, m,
                                                car[0], car[1], car[2]);

        if (!ce::approx(rec_car, result, tolerance))
        {
            cout << "cartesian table failed";
            cout << endl;
            checksum += 1;
        }
    }
    
    // Check functions against each other

    int const dimension = 3;
    int const number_of_tests = 100000;
    int const max_l = 8;
    
    Random_Number_Generator<int> rng_l(0, // lower
                                       max_l, // upper
                                       1492); // seed
    Random_Number_Generator<int> rng_m(-max_l,
                                       max_l,
                                       938);
    Random_Number_Generator<double> rng_x(-1,
                                          1, 
                                          582);
    
    for (int i = 0; i < number_of_tests; ++i)
    {
        int l = rng_l.scalar();
        int m = rng_l.scalar();
        
        while(std::abs(m) > l)
        {
            m = rng_m.scalar();
        }
        
        vector<double> cart = rng_x.vector(dimension);
        cart = vf::normalize(cart);
        vector<double> sph(2);
        mf::xyz_to_spherical(cart[0], cart[1], cart[2],
                             sph[0], sph[1]);
        
        // Check reverse back to cartesian
        
        vector<double> cart_new(dimension);
        mf::spherical_to_xyz(sph[0], sph[1],
                             cart_new[0], cart_new[1], cart_new[2]);
        
        if (!ce::approx(cart, cart_new, 1000 * tolerance))
        {
            cout << "cartesian reverse failed for case ";
            cout << i;
            cout << endl;
            cout << "\texpected:   ";
            cout << cart[0] << " " << cart[1] << " " << cart[2];
            cout << endl;
            cout << "\tcalculated: ";
            cout << cart_new[0] << " " << cart_new[1] << " " << cart_new[2];
            cout << endl;
            checksum += 1;
        }
        
        // Check spherical harmonic functions
        
        double rec_car = mf::spherical_harmonic_rec(l, m,
                                                    sph[0], sph[1]);
        double rec_sph = mf::spherical_harmonic_rec(l, m,
                                                    cart[0], cart[1], cart[2]);
        double tab_car = mf::spherical_harmonic(l, m,
                                                cart[0], cart[1], cart[2]);

        if (!(ce::approx(rec_car, rec_sph, 1000 * tolerance)
              && ce::approx(rec_car, tab_car, 1000 * tolerance)))
        {
            cout << "spherical harmonic comparison failed for case ";
            cout << i;
            cout << endl;
            cout << "\trec_car: ";
            cout << rec_car;
            cout << "\trec_sph: ";
            cout << rec_sph;
            cout << "\ttab_car: ";
            cout << tab_car;
            checksum += 1;
        }
    }
    
    return checksum;
}

int main()
{
    int checksum = 0;

    checksum += test_factorial();
    checksum += test_legendre_p();
    checksum += test_spherical_h();

    return checksum;
}

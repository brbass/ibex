#include <iostream>
#include <limits>

#include "Check_Equality.hh"
#include "Random_Number_Generator.hh"
#include "Vector_Functions.hh"
#include "Vector_Functions_2D.hh"
#include "Vector_Functions_3D.hh"

namespace ce = Check_Equality;
namespace vf = Vector_Functions;
namespace vf2 = Vector_Functions_2D;
namespace vf3 = Vector_Functions_3D;

using namespace std;

int test_vector_functions_2d()
{
    int checksum = 0;

    double tol = 10 * numeric_limits<double>::epsilon();
    
    vector<double> const t = {3, 1, 2, -7};
    vector<double> const x = {2, 4};
    vector<double> const y = {-3, 7};
    double const s = 2;

    vector<double> const add_xy = {-1, 11};
    if (!ce::approx(add_xy,
                    vf2::add(x, y),
                    tol))
    {
        cout << "vf2 add failed" << endl;
        checksum += 1;
    }
    if (!ce::approx(add_xy,
                    vf::add(x, y),
                    tol))
    {
        cout << "vf add failed" << endl;
        checksum += 1;
    }

    vector<double> const subtract_xy = {5, -3};
    if (!ce::approx(subtract_xy,
                    vf2::subtract(x, y),
                    tol))
    {
        cout << "vf2 subtract failed" << endl;
        checksum += 1;
    }
    if (!ce::approx(subtract_xy,
                    vf::subtract(x, y),
                    tol))
    {
        cout << "vf subtract failed" << endl;
        checksum += 1;
    }

    vector<double> const multiply_xs = {4, 8};
    if (!ce::approx(multiply_xs,
                    vf2::multiply(x, s),
                    tol))
    {
        cout << "vf2 multiply failed" << endl;
        checksum += 1;
    }
    if (!ce::approx(multiply_xs,
                    vf::multiply(x, s),
                    tol))
    {
        cout << "vf2 multiply failed" << endl;
        checksum += 1;
    }

    vector<double> const power_xs = {4, 16};
    // Not implemented for vf2
    if (!ce::approx(power_xs,
                    vf::power(x, s),
                    tol))
    {
        cout << "vf power failed" << endl;
        checksum += 1;
    }
    
    double const dot_xy = 22;
    if (!ce::approx(dot_xy,
                    vf2::dot(x, y),
                    tol))
    {
        cout << "vf2 dot failed" << endl;
        checksum += 1;
    }
    if (!ce::approx(dot_xy,
                    vf::dot(x, y),
                    tol))
    {
        cout << "vf dot failed" << endl;
        checksum += 1;
    }

    double cross_xy = 26;
    if (!ce::approx(cross_xy,
                    vf2::cross(x, y),
                    tol))
    {
        cout << "vf2 cross failed" << endl;
        checksum += 1;
    }
    // Not implemented for vf
    
    vector<double> const tensor_dot_tx = {10, -24};
    // Not implemented for vf2
    if (!ce::approx(tensor_dot_tx,
                    vf::tensor_dot(t, x),
                    tol))
    {
        cout << "vf tensor_dot failed" << endl;
        checksum += 1;
    }
    
    vector<double> const tensor_product_xy = {-6, 14, -12, 28};
    // Not implemented for vf2
    if (!ce::approx(tensor_product_xy,
                    vf::tensor_product(x, y),
                    tol))
    {
        cout << "vf tensor_dot failed" << endl;
        checksum += 1;
    }

    double const magnitude_x = 2 * sqrt(5);
    if (!ce::approx(magnitude_x,
                    vf2::magnitude(x),
                    tol))
    {
        cout << "vf2 magnitude failed" << endl;
        checksum += 1;
    }
    if (!ce::approx(magnitude_x,
                    vf::magnitude(x),
                    tol))
    {
        cout << "vf magnitude failed" << endl;
        checksum += 1;
    }

    double const magnitude_squared_x = 20;
    if (!ce::approx(magnitude_squared_x,
                    vf2::magnitude_squared(x),
                    tol))
    {
        cout << "vf2 magnitude_squared failed" << endl;
        checksum += 1;
    }
    if (!ce::approx(magnitude_squared_x,
                    vf::magnitude_squared(x),
                    tol))
    {
        cout << "vf magnitude_squared failed" << endl;
        checksum += 1;
    }
    
    vector<double> const normalize_x = {1. / sqrt(5.), 2. / sqrt(5.)};
    if (!ce::approx(normalize_x,
                    vf2::normalize(x),
                    tol))
    {
        cout << "vf2 normalization failed" << endl;
        checksum += 1;
    }
    if (!ce::approx(normalize_x,
                    vf::normalize(x),
                    tol))
    {
        cout << "vf normalization failed" << endl;
        checksum += 1;
    }
    
    return checksum;
}

int test_vector_functions_3d()
{
    int checksum = 0;

    double tol = 10 * numeric_limits<double>::epsilon();
    
    vector<double> const t = {-6, 6, 4, -1, -8, 0, 7, -10, -8};
    vector<double> const x = {-3, 4, -2};
    vector<double> const y = {9, -6, 7};
    double const s = -3;

    vector<double> const add_xy = {6, -2, 5};
    if (!ce::approx(add_xy,
                    vf3::add(x, y),
                    tol))
    {
        cout << "vf3 add failed" << endl;
        checksum += 1;
    }
    if (!ce::approx(add_xy,
                    vf::add(x, y),
                    tol))
    {
        cout << "vf add failed" << endl;
        checksum += 1;
    }

    vector<double> const subtract_xy = {-12, 10, -9};
    if (!ce::approx(subtract_xy,
                    vf3::subtract(x, y),
                    tol))
    {
        cout << "vf3 subtract failed" << endl;
        checksum += 1;
    }
    if (!ce::approx(subtract_xy,
                    vf::subtract(x, y),
                    tol))
    {
        cout << "vf subtract failed" << endl;
        checksum += 1;
    }

    vector<double> const multiply_xs = {9, -12, 6};
    if (!ce::approx(multiply_xs,
                    vf3::multiply(x, s),
                    tol))
    {
        cout << "vf3 multiply failed" << endl;
        checksum += 1;
    }
    if (!ce::approx(multiply_xs,
                    vf::multiply(x, s),
                    tol))
    {
        cout << "vf3 multiply failed" << endl;
        checksum += 1;
    }

    vector<double> const power_xs = {-1. / 27., 1. / 64., -1. / 8.};
    // Not implemented for vf3
    if (!ce::approx(power_xs,
                    vf::power(x, s),
                    tol))
    {
        cout << "vf power failed" << endl;
        checksum += 1;
    }
    
    double const dot_xy = -65;
    if (!ce::approx(dot_xy,
                    vf3::dot(x, y),
                    tol))
    {
        cout << "vf3 dot failed" << endl;
        checksum += 1;
    }
    if (!ce::approx(dot_xy,
                    vf::dot(x, y),
                    tol))
    {
        cout << "vf dot failed" << endl;
        checksum += 1;
    }

    vector<double> const cross_xy = {34, -29, -45};
    if (!ce::approx(cross_xy,
                    vf3::cross(x, y),
                    tol))
    {
        cout << "vf3 cross failed" << endl;
        checksum += 1;
    }
    // Not implemented for vf
    
    vector<double> const tensor_dot_tx = {16, 3, -18};
    // Not implemented for vf3
    if (!ce::approx(tensor_dot_tx,
                    vf::tensor_dot(t, x),
                    tol))
    {
        cout << "vf tensor_dot failed" << endl;
        checksum += 1;
    }
    
    vector<double> const tensor_product_xy = {-27, 18, -21, 36, -24, 28, -18, 12, -14};
    // Not implemented for vf3
    if (!ce::approx(tensor_product_xy,
                    vf::tensor_product(x, y),
                    tol))
    {
        cout << "vf tensor_dot failed" << endl;
        checksum += 1;
    }

    double const magnitude_x = sqrt(29);
    if (!ce::approx(magnitude_x,
                    vf3::magnitude(x),
                    tol))
    {
        cout << "vf3 magnitude failed" << endl;
        checksum += 1;
    }
    if (!ce::approx(magnitude_x,
                    vf::magnitude(x),
                    tol))
    {
        cout << "vf magnitude failed" << endl;
        checksum += 1;
    }

    double const magnitude_squared_x = 29;
    if (!ce::approx(magnitude_squared_x,
                    vf3::magnitude_squared(x),
                    tol))
    {
        cout << "vf3 magnitude_squared failed" << endl;
        checksum += 1;
    }
    if (!ce::approx(magnitude_squared_x,
                    vf::magnitude_squared(x),
                    tol))
    {
        cout << "vf magnitude_squared failed" << endl;
        checksum += 1;
    }
    
    vector<double> const normalize_x = {-3. / sqrt(29), 4. / sqrt(29), -2. / sqrt(29)};
    if (!ce::approx(normalize_x,
                    vf3::normalize(x),
                    tol))
    {
        cout << "vf3 normalization failed" << endl;
        checksum += 1;
    }
    if (!ce::approx(normalize_x,
                    vf::normalize(x),
                    tol))
    {
        cout << "vf normalization failed" << endl;
        checksum += 1;
    }
    
    return checksum;
}

int main()
{
    int checksum = 0;

    checksum += test_vector_functions_2d();
    checksum += test_vector_functions_3d();
}

#include "LDFE_Quadrature.hh"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>

#include "Check.hh"
#include "XML_Functions.hh"

using namespace std;

namespace // anonymous
{
    int get_number_of_ordinates(int dimension,
                                int rule)
    {
        vector<int> rules = {4, 16, 64};
    
        switch(rule)
        {
        case 1:
            return rules[0] * 4 * (dimension - 1);
        case 2:
            return rules[1] * 4 * (dimension - 1);
        case 3:
            return rules[2] * 4 * (dimension - 1);
        default:
            AssertMsg(false, "rule not found");
            return -1;
        }
    }
}


LDFE_Quadrature::
LDFE_Quadrature(int dimension,
                int number_of_scattering_moments,
                int rule):
    Angular_Discretization(dimension,
                           number_of_scattering_moments,
                           get_number_of_ordinates(dimension,
                                                   rule)),
    rule_(rule)
{
    reflection_tolerance_ = 1000 * numeric_limits<double>::epsilon();
    initialize_quadrature();

    directions_.resize(number_of_ordinates_);
    for (int o = 0; o < number_of_ordinates_; ++o)
    {
        vector<double> direction(dimension_);
        
        for (int d = 0; d < dimension_; ++d)
        {
            direction[d] = ordinates_[d + dimension * o];
        }

        directions_[o] = direction;
    }
    
    check_class_invariants();
}

void LDFE_Quadrature::
initialize_quadrature()
{
    switch(rule_)
    {
    case 1:
        initialize_1();
        break;
    case 2:
        initialize_2();
        break;
    case 3: 
        initialize_3();
        break;
    }
    
    switch(dimension_)
    {
    case 2:
        for (int i = 0; i < number_of_ordinates_; ++i)
        {
            ordinates_.push_back(mu_[i]);
            ordinates_.push_back(eta_[i]);
        }
        break;
    case 3:
        for (int i = 0; i < number_of_ordinates_; ++i)
        {
            ordinates_.push_back(mu_[i]);
            ordinates_.push_back(eta_[i]);
            ordinates_.push_back(xi_[i]);
        }
        break;
    }
}

void LDFE_Quadrature::
check_class_invariants() const
{
    Assert(mu_.size() == number_of_ordinates_);
    Assert(eta_.size() == number_of_ordinates_);
    Assert(xi_.size() == number_of_ordinates_);
    Assert(weights_.size() == number_of_ordinates_);
    Assert(ordinates_.size() == number_of_ordinates_ * dimension_);
    Assert(directions_.size() == number_of_ordinates_);
}

void LDFE_Quadrature::
output(pugi::xml_node &output_node) const
{
    pugi::xml_node ldfe = output_node.append_child("angular_discretization");

    XML_Functions::append_child(ldfe, "ldfe_quadrature", "quadrature_type");
    XML_Functions::append_child(ldfe, dimension_, "dimension");
    XML_Functions::append_child(ldfe, rule_, "rule");
    XML_Functions::append_child(ldfe, number_of_moments_, "number_of_moments");
    XML_Functions::append_child(ldfe, number_of_ordinates_, "number_of_ordinates");
    XML_Functions::append_child(ldfe, ordinates_, "ordinates", "dimension-ordinate");
    XML_Functions::append_child(ldfe, weights_, "weights", "ordinate");
    XML_Functions::append_child(ldfe, mu_, "mu", "ordinate");
    XML_Functions::append_child(ldfe, eta_, "eta", "ordinate");
    XML_Functions::append_child(ldfe, xi_, "xi", "ordinate");
    XML_Functions::append_child(ldfe, weights_, "weights", "ordinate");
}

int LDFE_Quadrature::
reflect_ordinate(int o,
                 vector<double> const &normal) const
{
    switch(dimension_)
    {
    case 2:
        if (abs(abs(normal[0]) - 1) < reflection_tolerance_)
        {
            if (mu_[o] > 0)
            {
                return o - 2;
            }
            else
            {
                return o + 2;
            }
        }
        else if (abs(abs(normal[1]) - 1) < reflection_tolerance_)
        {
            if (eta_[o] > 0)
            {
                return o - 1;
            }
            else
            {
                return o + 1;
            }
        }
        else 
        {
            AssertMsg(false, "surface must be in x, y, or z plane: {" + to_string(normal[0]) + ", " + to_string(normal[1]) + "}");
            
            return o;
        }
    case 3:
        if (abs(abs(normal[0]) - 1) < reflection_tolerance_)
        {
            if (mu_[o] > 0)
            {
                return o - 4;
            }
            else
            {
                return o + 4;
            }
        }
        else if (abs(abs(normal[1]) - 1) < reflection_tolerance_)
        {
            if (eta_[o] > 0)
            {
                return o - 2;
            }
            else
            {
                return o + 2;
            }
        }
        else if (abs(abs(normal[2]) - 1) < reflection_tolerance_)
        {
            if (xi_[o] > 0)
            {
                return o - 1;
            }
            else
            {
                return o + 1;
            }
        }
        else 
        {
            AssertMsg(false, "surface must be in x, y, or z plane: {" + to_string(normal[0]) + ", " + to_string(normal[1]) + ", " + to_string(normal[2]) + "}");
            
            return o;
        }
    default:
        AssertMsg(false, "dimension not found");

        return o;
    }
}

void LDFE_Quadrature::
initialize_1()
{
    int num_vals = 4;
    int number_of_ordinates_per_octant = 4;
    
    vector<double> vals
        = {2.09769106510711E-01, 2.09769106510711E-01, 9.54983687770317E-01, 3.39836909454130E-01,
           9.54983687770317E-01, 2.09769106510711E-01, 2.09769106510711E-01, 3.39836909454103E-01,
           2.09769106510711E-01, 9.54983687770317E-01, 2.09769106510711E-01, 3.39836909454103E-01,
           5.77350269189625E-01, 5.77350269189625E-01, 5.77350269189625E-01, 5.51285598432481E-01};

    mu_.resize(0);
    eta_.resize(0);
    xi_.resize(0);
    
    switch(dimension_)
    {
    case 2:
        for (int i = 0; i < number_of_ordinates_per_octant; ++i)
        {
            for (int mu_mult = -1; mu_mult < 2; mu_mult += 2)
            {
                for (int eta_mult = -1; eta_mult < 2; eta_mult += 2)
                {
                    int xi_mult = 1;
                    
                    mu_.push_back(vals[0 + num_vals * i] * mu_mult);
                    eta_.push_back(vals[1 + num_vals * i] * eta_mult);
                    xi_.push_back(vals[2 + num_vals * i] * xi_mult);
                    weights_.push_back(vals[3 + num_vals * i]);
                }
            }
        }
        break;
    case 3:
        for (int i = 0; i < number_of_ordinates_per_octant; ++i)
        {
            for (int mu_mult = -1; mu_mult < 2; mu_mult += 2)
            {
                for (int eta_mult = -1; eta_mult < 2; eta_mult += 2)
                {
                    for (int xi_mult = -1; xi_mult < 2; xi_mult += 2)
                    {
                        mu_.push_back(vals[0 + num_vals * i] * mu_mult);
                        eta_.push_back(vals[1 + num_vals * i] * eta_mult);
                        xi_.push_back(vals[2 + num_vals * i] * xi_mult);
                        weights_.push_back(vals[3 + num_vals * i]);
                    }
                }
            }
        }
        break;
    }
}

void LDFE_Quadrature::
initialize_2()
{
    int num_vals = 4;
    int number_of_ordinates_per_octant = 16;

    vector<double> vals
        = {1.22785291159728E-01, 1.22785291159728E-01, 9.84808379609780E-01, 5.26559082615718E-02,
           5.31419325021509E-01, 1.03129501638898E-01, 8.40807829938206E-01, 9.95720041972116E-02,
           1.03129501638901E-01, 5.31419325021504E-01, 8.40807829938209E-01, 9.95720041972139E-02,
           2.35702260395515E-01, 2.35702260395515E-01, 9.42809041582063E-01, 8.80369927980968E-02,
           8.40807829938193E-01, 1.03129501638887E-01, 5.31419325021532E-01, 9.95720041972128E-02,
           9.84808379609774E-01, 1.22785291159753E-01, 1.22785291159753E-01, 5.26559082615709E-02,
           8.40807829938210E-01, 5.31419325021502E-01, 1.03129501638902E-01, 9.95720041972133E-02,
           9.42809041582063E-01, 2.35702260395515E-01, 2.35702260395515E-01, 8.80369927981024E-02,
           1.03129501638887E-01, 8.40807829938193E-01, 5.31419325021531E-01, 9.95720041972128E-02,
           5.31419325021508E-01, 8.40807829938206E-01, 1.03129501638899E-01, 9.95720041972108E-02,
           1.22785291159756E-01, 9.84808379609773E-01, 1.22785291159756E-01, 5.26559082615713E-02,
           2.35702260395515E-01, 9.42809041582063E-01, 2.35702260395515E-01, 8.80369927981055E-02,
           6.86947072008541E-01, 6.86947072008541E-01, 2.37081084268193E-01, 1.32024927825519E-01,
           2.37081084268253E-01, 6.86947072008531E-01, 6.86947072008531E-01, 1.32024927825525E-01,
           6.86947072008531E-01, 2.37081084268254E-01, 6.86947072008531E-01, 1.32024927825526E-01,
           5.77350269189625E-01, 5.77350269189625E-01, 5.77350269189625E-01, 1.55210814955923E-01};
    
    mu_.resize(0);
    eta_.resize(0);
    xi_.resize(0);
    
    switch(dimension_)
    {
    case 2:
        for (int i = 0; i < number_of_ordinates_per_octant; ++i)
        {
            for (int mu_mult = -1; mu_mult < 2; mu_mult += 2)
            {
                for (int eta_mult = -1; eta_mult < 2; eta_mult += 2)
                {
                    int xi_mult = 1;
                    
                    mu_.push_back(vals[0 + num_vals * i] * mu_mult);
                    eta_.push_back(vals[1 + num_vals * i] * eta_mult);
                    xi_.push_back(vals[2 + num_vals * i] * xi_mult);
                    weights_.push_back(vals[3 + num_vals * i]);
                }
            }
        }
        break;
    case 3:
        for (int i = 0; i < number_of_ordinates_per_octant; ++i)
        {
            for (int mu_mult = -1; mu_mult < 2; mu_mult += 2)
            {
                for (int eta_mult = -1; eta_mult < 2; eta_mult += 2)
                {
                    for (int xi_mult = -1; xi_mult < 2; xi_mult += 2)
                    {
                        mu_.push_back(vals[0 + num_vals * i] * mu_mult);
                        eta_.push_back(vals[1 + num_vals * i] * eta_mult);
                        xi_.push_back(vals[2 + num_vals * i] * xi_mult);
                        weights_.push_back(vals[3 + num_vals * i]);
                    }
                }
            }
        }
    }
}

void LDFE_Quadrature::
initialize_3()
{
    int num_vals = 4;
    int number_of_ordinates_per_octant = 64;

    vector<double> vals =
        {4.38752387852336E-02, 4.38752387852336E-02, 9.98073106963150E-01, 1.01011818772561E-02,
         2.23194054503770E-01, 4.34603636088626E-02, 9.73804708773352E-01, 1.47114091091454E-02,
         4.34603636088619E-02, 2.23194054503772E-01, 9.73804708773351E-01, 1.47114091091449E-02,
         9.90147542976669E-02, 9.90147542976669E-02, 9.90147542976674E-01, 1.31319081660192E-02,
         3.91754724400026E-01, 5.13254348711689E-02, 9.18636998844235E-01, 2.01635596440037E-02,
         6.25662723829225E-01, 5.09223830015233E-02, 7.78429872833796E-01, 2.42900177326475E-02,
         4.50060753113050E-01, 2.86782825010254E-01, 8.45695530191836E-01, 3.00096725803229E-02,
         4.92365963917330E-01, 1.23091490979332E-01, 8.61640436855329E-01, 2.51087542402306E-02,
         5.13254348711655E-02, 3.91754724400022E-01, 9.18636998844237E-01, 2.01635596440024E-02,
         2.86782825010250E-01, 4.50060753113052E-01, 8.45695530191837E-01, 3.00096725803234E-02,
         5.09223830015209E-02, 6.25662723829230E-01, 7.78429872833792E-01, 2.42900177326465E-02,
         1.23091490979332E-01, 4.92365963917330E-01, 8.61640436855329E-01, 2.51087542402332E-02,
         3.35052576156500E-01, 3.35052576156500E-01, 8.80613162757510E-01, 2.78797764691254E-02,
         9.47958080258904E-02, 2.86222206910509E-01, 9.53462428757418E-01, 1.90475751591700E-02,
         2.86222206910502E-01, 9.47958080259107E-02, 9.53462428757418E-01, 1.90475751591731E-02,
         2.35702260395515E-01, 2.35702260395515E-01, 9.42809041582063E-01, 2.20620660106520E-02,
         7.78429872833780E-01, 5.09223830015111E-02, 6.25662723829246E-01, 2.42900177326465E-02,
         9.18636998844232E-01, 5.13254348711738E-02, 3.91754724400033E-01, 2.01635596440023E-02,
         8.45695530191835E-01, 2.86782825010258E-01, 4.50060753113049E-01, 3.00096725803224E-02,
         8.61640436855329E-01, 1.23091490979332E-01, 4.92365963917330E-01, 2.51087542402348E-02,
         9.73804708773349E-01, 4.34603636088579E-02, 2.23194054503781E-01, 1.47114091091450E-02,
         9.98073106963150E-01, 4.38752387852414E-02, 4.38752387852413E-02, 1.01011818772558E-02,
         9.73804708773351E-01, 2.23194054503774E-01, 4.34603636088611E-02, 1.47114091091441E-02,
         9.90147542976674E-01, 9.90147542976673E-02, 9.90147542976674E-02, 1.31319081660202E-02,
         8.45695530191832E-01, 4.50060753113044E-01, 2.86782825010275E-01, 3.00096725803228E-02,
         9.18636998844233E-01, 3.91754724400032E-01, 5.13254348711725E-02, 2.01635596440023E-02,
         7.78429872833787E-01, 6.25662723829236E-01, 5.09223830015169E-02, 2.42900177326455E-02,
         8.61640436855329E-01, 4.92365963917330E-01, 1.23091490979332E-01, 2.51087542402358E-02,
         9.53462428757418E-01, 2.86222206910519E-01, 9.47958080258602E-02, 1.90475751591677E-02,
         8.80613162757525E-01, 3.35052576156479E-01, 3.35052576156479E-01, 2.78797764691281E-02,
         9.53462428757418E-01, 9.47958080258961E-02, 2.86222206910507E-01, 1.90475751591710E-02,
         9.42809041582063E-01, 2.35702260395515E-01, 2.35702260395515E-01, 2.20620660106529E-02,
         5.09223830015081E-02, 7.78429872833776E-01, 6.25662723829251E-01, 2.42900177326452E-02,
         2.86782825010257E-01, 8.45695530191835E-01, 4.50060753113050E-01, 3.00096725803224E-02,
         5.13254348711753E-02, 9.18636998844231E-01, 3.91754724400035E-01, 2.01635596440032E-02,
         1.23091490979332E-01, 8.61640436855329E-01, 4.92365963917330E-01, 2.51087542402352E-02,
         4.50060753113045E-01, 8.45695530191833E-01, 2.86782825010272E-01, 3.00096725803231E-02,
         6.25662723829232E-01, 7.78429872833790E-01, 5.09223830015193E-02, 2.42900177326471E-02,
         3.91754724400031E-01, 9.18636998844233E-01, 5.13254348711724E-02, 2.01635596440025E-02,
         4.92365963917331E-01, 8.61640436855328E-01, 1.23091490979332E-01, 2.51087542402323E-02,
         4.34603636088583E-02, 9.73804708773350E-01, 2.23194054503780E-01, 1.47114091091460E-02,
         2.23194054503777E-01, 9.73804708773350E-01, 4.34603636088597E-02, 1.47114091091451E-02,
         4.38752387852359E-02, 9.98073106963150E-01, 4.38752387852358E-02, 1.01011818772555E-02,
         9.90147542976673E-02, 9.90147542976674E-01, 9.90147542976674E-02, 1.31319081660197E-02,
         2.86222206910522E-01, 9.53462428757418E-01, 9.47958080258519E-02, 1.90475751591660E-02,
         9.47958080258752E-02, 9.53462428757418E-01, 2.86222206910514E-01, 1.90475751591669E-02,
         3.35052576156480E-01, 8.80613162757525E-01, 3.35052576156480E-01, 2.78797764691266E-02,
         2.35702260395515E-01, 9.42809041582063E-01, 2.35702260395515E-01, 2.20620660106607E-02,
         7.02505464432133E-01, 7.02505464432132E-01, 1.13895324249885E-01, 2.78705916500711E-02,
         5.29915360873601E-01, 7.69766972496521E-01, 3.55877111323192E-01, 3.51928031147286E-02,
         7.69766972496521E-01, 5.29915360873602E-01, 3.55877111323192E-01, 3.51928031147287E-02,
         6.80413817439771E-01, 6.80413817439771E-01, 2.72165526975908E-01, 3.37687299459803E-02,
         3.55877111323199E-01, 7.69766972496527E-01, 5.29915360873588E-01, 3.51928031147297E-02,
         1.13895324249909E-01, 7.02505464432130E-01, 7.02505464432130E-01, 2.78705916500730E-02,
         3.55877111323192E-01, 5.29915360873602E-01, 7.69766972496521E-01, 3.51928031147292E-02,
         2.72165526975908E-01, 6.80413817439771E-01, 6.80413817439771E-01, 3.37687299459794E-02,
         7.69766972496528E-01, 3.55877111323199E-01, 5.29915360873587E-01, 3.51928031147294E-02,
         5.29915360873602E-01, 3.55877111323191E-01, 7.69766972496521E-01, 3.51928031147287E-02,
         7.02505464432131E-01, 1.13895324249906E-01, 7.02505464432131E-01, 2.78705916500716E-02,
         6.80413817439771E-01, 2.72165526975908E-01, 6.80413817439771E-01, 3.37687299459807E-02,
         4.85081580043393E-01, 4.85081580043393E-01, 7.27593101537672E-01, 3.83614334577098E-02,
         7.27593101537650E-01, 4.85081580043409E-01, 4.85081580043409E-01, 3.83614334577117E-02,
         4.85081580043406E-01, 7.27593101537654E-01, 4.85081580043406E-01, 3.83614334577106E-02,
         5.77350269189625E-01, 5.77350269189625E-01, 5.77350269189625E-01, 4.01265145828269E-02};

    mu_.resize(0);
    eta_.resize(0);
    xi_.resize(0);


    
    switch(dimension_)
    {
    case 2:
        for (int i = 0; i < number_of_ordinates_per_octant; ++i)
        {
            for (int mu_mult = -1; mu_mult < 2; mu_mult += 2)
            {
                for (int eta_mult = -1; eta_mult < 2; eta_mult += 2)
                {
                    int xi_mult = 1;
                    
                    mu_.push_back(vals[0 + num_vals * i] * mu_mult);
                    eta_.push_back(vals[1 + num_vals * i] * eta_mult);
                    xi_.push_back(vals[2 + num_vals * i] * xi_mult);
                    weights_.push_back(vals[3 + num_vals * i]);
                }
            }
        }
        break;
    case 3:
        for (int i = 0; i < number_of_ordinates_per_octant; ++i)
        {
            for (int mu_mult = -1; mu_mult < 2; mu_mult += 2)
            {
                for (int eta_mult = -1; eta_mult < 2; eta_mult += 2)
                {
                    for (int xi_mult = -1; xi_mult < 2; xi_mult += 2)
                    {
                        mu_.push_back(vals[0 + num_vals * i] * mu_mult);
                        eta_.push_back(vals[1 + num_vals * i] * eta_mult);
                        xi_.push_back(vals[2 + num_vals * i] * xi_mult);
                        weights_.push_back(vals[3 + num_vals * i]);
                    }
                }
            }
        }
    }
}


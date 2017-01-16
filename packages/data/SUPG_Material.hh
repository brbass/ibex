#ifndef SUPG_Material_hh
#define SUPG_Material_hh

#include "Material.hh"

/*
  Holds cross sections with dimensionally-dependent moments
*/
class SUPG_Material : public Material
{
public:

    SUPG_Material(int index,
                  std::shared_ptr<Angular_Discretization> angular_discretization,
                  std::shared_ptr<Energy_Discretization> energy_discretization,
                  std::vector<double> const &sigma_t,
                  std::vector<double> const &sigma_s,
                  std::vector<double> const &nu,
                  std::vector<double> const &sigma_f,
                  std::vector<double> const &chi,
                  std::vector<double> const &internal_source);

    // Check class invariants
    virtual void check_class_invariants() const override;

    // Output data to XML file
    virtual void output(XML_Node output_node) const override;
};

#endif

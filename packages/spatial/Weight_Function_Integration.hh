#ifndef Weight_Function_Integration_hh
#define Weight_Function_Integration_hh

#include <memory>
#include <vector>

#include "Integration_Mesh.hh"
#include "Weight_Function.hh"

class Angular_Discretization;
class Basis_Function;
class Energy_Discretization;
class Material;
class Weak_Spatial_Discretization_Options;

class Weight_Function_Integration
{
public:

    struct Material_Data
    {
        std::vector<double> sigma_t;
        std::vector<double> sigma_s;
        std::vector<double> nu;
        std::vector<double> sigma_f;
        std::vector<double> chi;
        std::vector<double> internal_source;
        std::vector<double> norm;
        std::vector<double> boundary_sources;
    };
    
    Weight_Function_Integration(int number_of_points,
                                std::shared_ptr<Weak_Spatial_Discretization_Options> weak_options,
                                std::vector<std::shared_ptr<Basis_Function> > const &bases,
                                std::vector<std::shared_ptr<Weight_Function> > const &weights);
    
    // Perform integration and put result into weight functions
    void perform_integration();
    
private:

    // Put volume, surface and material integrals into weight functions
    void put_integrals_into_weight(std::vector<Weight_Function::Integrals> const &integrals,
                                   std::vector<Material_Data> const &materials);
    
    // Perform all volume integrals
    void perform_volume_integration(std::vector<Weight_Function::Integrals> &integrals,
                                    std::vector<Material_Data> &materials) const;

    // Normalize the material integrals, if applicable
    void normalize_materials(std::vector<Material_Data> &materials) const;
    
    // Add weight function cell values to global integrals
    void add_volume_weight(std::shared_ptr<Integration_Mesh::Cell> const cell,
                           double quad_weight,
                           std::vector<double> const &w_val,
                           std::vector<std::vector<double> > const &w_grad,
                           std::vector<Weight_Function::Integrals> &integrals) const;

    // Add basis/weight function cell values to global integrals
    void add_volume_basis_weight(std::shared_ptr<Integration_Mesh::Cell> const cell,
                                 double quad_weight,
                                 std::vector<double> const &b_val,
                                 std::vector<std::vector<double> > const &b_grad,
                                 std::vector<double> const &w_val,
                                 std::vector<std::vector<double> > const &w_grad,
                                 std::vector<std::vector<int> > const &weight_basis_indices,
                                 std::vector<Weight_Function::Integrals> &integrals) const;

    // Add material cell values to global integrals
    void add_volume_material(std::shared_ptr<Integration_Mesh::Cell> const cell,
                             double quad_weight,
                             std::vector<double> const &b_val,
                             std::vector<double> const &w_val,
                             std::vector<std::vector<double> > const &w_grad,
                             std::vector<std::vector<int> > const &weight_basis_indices,
                             std::shared_ptr<Material> point_material,
                             std::vector<Material_Data> &materials) const;

    // Add local vector to total vector
    void add_vector_to_total(std::vector<double> const &local_data,
                             std::vector<double> &data);
    
    // Add local integrals to global integrals (for parallel)
    void add_integral_sets(std::vector<Weight_Function::Integrals> const &local_integrals,
                           std::vector<Weight_Function::Integrals> &integrals);

    // Add local materials to global materials (for parallel)
    void add_material_sets(std::vector<Material_Data> const &local_materials,
                           std::vector<Material_Data> &materials);
    
    // Perform all surface integrals
    void perform_surface_integration(std::vector<Weight_Function::Integrals> &integrals,
                                     std::vector<Material_Data> &materials) const;

    // Add boundary source integral to global integrals
    void add_surface_source(std::shared_ptr<Integration_Mesh::Surface> const surface,
                            double quad_weight,
                            std::vector<double> const &w_val,
                            std::vector<int> const &weight_surface_indices,
                            std::shared_ptr<Boundary_Source> boundary_source,
                            std::vector<Material_Data> &materials) const;
    
    // Add weight function surface values to global integrals
    void add_surface_weight(std::shared_ptr<Integration_Mesh::Surface> const surface,
                            double quad_weight,
                            std::vector<double> const &w_val,
                            std::vector<int> const &weight_surface_indices,
                            std::vector<Weight_Function::Integrals> &integrals) const;
    
    // Add basis/weight function surface values to global integrals
    void add_surface_basis_weight(std::shared_ptr<Integration_Mesh::Surface> const surface,
                                  double quad_weight,
                                  std::vector<double> const &b_val,
                                  std::vector<double> const &w_val,
                                  std::vector<int> const &weight_surface_indices,
                                  std::vector<std::vector<int> > const &weight_basis_indices,
                                  std::vector<Weight_Function::Integrals> &integrals) const;
    
    // Get cross section data
    void get_cross_sections(std::shared_ptr<Material> material,
                            std::vector<double> &sigma_t,
                            std::vector<double> &sigma_s,
                            std::vector<double> &chi_nu_sigma_f,
                            std::vector<double> &internal_source) const;
    
    // Get a material from the material data
    void get_material(int index,
                      Material_Data const &material_data,
                      std::shared_ptr<Material> &material) const;

    // Get boundary sources from the material data
    void get_boundary_sources(int index,
                              Material_Data const &material_data,
                              std::vector<std::shared_ptr<Boundary_Source> > &boundary_sources) const;
    
    // Get the flux for a specific point, given basis coefficients
    void get_flux(std::shared_ptr<Integration_Mesh::Cell> const cell,
                  std::vector<double> const &b_val,
                  std::vector<double> &flux) const;
    
    // Initialize material data to zero
    void initialize_materials(std::vector<Material_Data> &materials) const;

    // Initialize integral data to zero
    void initialize_integrals(std::vector<Weight_Function::Integrals> &integrals) const;

    // Data
    int dimension_;
    int number_of_points_;
    std::shared_ptr<Weak_Spatial_Discretization_Options> options_;
    std::vector<std::shared_ptr<Basis_Function> > bases_;
    std::vector<std::shared_ptr<Weight_Function> > weights_;
    std::shared_ptr<Solid_Geometry> solid_;
    std::shared_ptr<Angular_Discretization> angular_;
    std::shared_ptr<Energy_Discretization> energy_;
    std::shared_ptr<Integration_Mesh> mesh_;
};
    
#endif

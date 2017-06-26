from mesh_functions import *

# Get points and find connectivity for pincell
def get_slab_discretization(length,
                            num_points,
                            num_neighbors_basis,
                            num_neighbors_weight,
                            use_constant_radius):
    points = get_cartesian_points_1d(length,
                                     num_points)
    
    return get_connectivity(num_points,
                            points,
                            num_neighbors_basis,
                            num_neighbors_weight,
                            use_constant_radius)

def output_slab_discretization(length,
                               num_points,
                               num_neighbors_basis,
                               num_neighbors_weight,
                               use_constant_radius):
    num_points, points, radius_basis, radius_weight, b_neighbors, b_neighbor_distances, w_neighbors, w_neighbor_distances, bw_neighbors, bw_neighbor_distances \
        = get_slab_discretization(length,
                                  num_points,
                                  num_neighbors_basis,
                                  num_neighbors_weight,
                                  use_constant_radius)
    output_path = "in/slab_{}_{}_{}_{}_{}.xml".format(length, num_points, num_neighbors_basis, num_neighbors_weight, int(use_constant_radius))

    # Get spatial discretization node
    node = xml_discretization(output_path,
                              num_points,
                              points,
                              radius_basis,
                              radius_weight,
                              b_neighbors,
                              b_neighbor_distances,
                              w_neighbors,
                              w_neighbor_distances,
                              bw_neighbors,
                              bw_neighbor_distances)

    # Add solid geometry information
    solid = et.SubElement(node, "solid_geometry")
    et.SubElement(solid, "dimension").text = "1"

    # Add surfaces
    surfaces = et.SubElement(solid, "surfaces")
    l2 = length / 2
    origins = np.array([-l2, l2])
    normals = np.array([-1, 1])
    boundary_sources = [0, 0]
    et.SubElement(surfaces, "number_of_surfaces").text = "2"
    for i in range(2):
        surface = et.SubElement(surfaces, "surface")
        surface.set("index", str(i))
        surface.set("shape", "plane")
        surface.set("type", "boundary")
        et.SubElement(surface, "origin").text = numpy_to_text(origins[i])
        et.SubElement(surface, "normal").text = numpy_to_text(normals[i])
        et.SubElement(surface, "boundary_source").text = str(boundary_sources[i])
        
    # Add regions
    regions = et.SubElement(solid, "regions")
    region = et.SubElement(regions, "region")
    region.set("index", "0")
    region.set("material", str(0))
    et.SubElement(regions, "number_of_regions").text = "1"
    for i in range(2):
        sr = et.SubElement(region, "surface_relation")
        sr.set("surface", str(i))
        sr.set("relation", "inside")
        
    et.ElementTree(node).write(output_path,
                               pretty_print=True,
                               xml_declaration=True)
    
if __name__ == '__main__':
    if (len(sys.argv) != 5):
        print("slab_functions.py [length num_points num_neighbors_basis num_neighbors_weight]")
        sys.exit()
    length = float(sys.argv[1])
    num_points = int(sys.argv[2])
    num_neighbors_basis = int(sys.argv[3])
    num_neighbors_weight = int(sys.argv[4])
    
    output_slab_discretization(length,
                               num_points,
                               num_neighbors_basis,
                               num_neighbors_weight,
                               True)
    
    

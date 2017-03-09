from mesh_functions import *

# Get points and find connectivity for pincell
def get_square_discretization(length,
                              num_points_xy,
                              num_neighbors_basis,
                              num_neighbors_weight):

    num_points, points = get_cartesian_points(length,
                                              num_points_xy)
    
    return get_connectivity(num_points,
                            points,
                            num_neighbors_basis,
                            num_neighbors_weight)

def output_square_discretization(length,
                                 num_points_xy,
                                 num_neighbors_basis,
                                 num_neighbors_weight):
    num_points, points, radius_basis, radius_weight, b_neighbors, b_neighbor_distances, w_neighbors, w_neighbor_distances, bw_neighbors, bw_neighbor_distances \
        = get_square_discretization(length,
                                    num_points_xy,
                                    num_neighbors_basis,
                                    num_neighbors_weight)
    
    output_path = "in/square_{}_{}_{}_{}.xml".format(length, num_points_xy, num_neighbors_basis, num_neighbors_weight)

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
    et.SubElement(solid, "dimension").text = "2"
    surfaces = et.SubElement(solid, "surfaces")
    l2 = length / 2
    origins = np.array([[-l2, 0.],
                        [l2, 0.],
                        [0., -l2],
                        [0., l2]])
    normals = np.array([[-1., 0.],
                        [1., 0.],
                        [0., -1],
                        [0., 1.]])
    boundary_sources = [0, 0, 0, 0]
    et.SubElement(surfaces, "number_of_surfaces").text = "4"
    for i in range(4):
        surface = et.SubElement(surfaces, "surface")
        surface.set("index", str(i))
        surface.set("shape", "plane")
        surface.set("type", "boundary")
        et.SubElement(surface, "origin").text = numpy_to_text(origins[i])
        et.SubElement(surface, "normal").text = numpy_to_text(normals[i])
        et.SubElement(surface, "boundary_source").text = str(boundary_sources[i])
        
    regions = et.SubElement(solid, "regions")
    region = et.SubElement(regions, "region")
    region.set("index", "0")
    region.set("material", str(0))
    et.SubElement(regions, "number_of_regions").text = "1"
    for i in range(4):
        sr = et.SubElement(region, "surface_relation")
        sr.set("surface", str(i))
        sr.set("relation", "inside")
        
    et.ElementTree(node).write(output_path,
                               pretty_print=True,
                               xml_declaration=True)

    
if __name__ == '__main__':
    if (len(sys.argv) != 5):
        print("square_functions.py [length num_points_xy num_neighbors_basis num_neighbors_weight]")
        sys.exit()
    length = float(sys.argv[1])
    num_points_xy = int(sys.argv[2])
    num_neighbors_basis = int(sys.argv[3])
    num_neighbors_weight = int(sys.argv[4])
    
    output_square_discretization(length,
                                 num_points_xy,
                                 num_neighbors_basis,
                                 num_neighbors_weight)
    

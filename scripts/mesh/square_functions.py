from mesh_functions import *
import argparse

# Get points and find connectivity for pincell
def get_square_discretization(length,
                              num_points_xy,
                              num_neighbors_basis,
                              num_neighbors_weight,
                              use_constant_radius):

    num_points, points = get_square(length,
                                    num_points_xy)
    
    return get_connectivity(num_points,
                            points,
                            num_neighbors_basis,
                            num_neighbors_weight,
                            use_constant_radius)

def output_square_discretization_points(length,
                                        num_points_xy,
                                        exclude_corners = False):
    output_path = "square_{}_{}.xml".format(length, num_points_xy)
    num_points, points = get_square(length,
                                    num_points_xy,
                                    exclude_corners)
    node = xml_points(2, # dimension
                      num_points,
                      points)
    spatial_node = node.find("spatial_discretization")
    
    et.ElementTree(node).write(output_path,
                               pretty_print=True,
                               xml_declaration=True)
    
    
    
def output_square_discretization(length,
                                 num_points_xy,
                                 num_neighbors_basis,
                                 num_neighbors_weight,
                                 use_constant_radius):
    num_points, points, radius_basis, radius_weight, b_neighbors, b_neighbor_distances, w_neighbors, w_neighbor_distances, bw_neighbors, bw_neighbor_distances \
        = get_square_discretization(length,
                                    num_points_xy,
                                    num_neighbors_basis,
                                    num_neighbors_weight,
                                    use_constant_radius)
    
    output_path = "square_{}_{}_{}_{}_{}.xml".format(length, num_points_xy, num_neighbors_basis, num_neighbors_weight, int(use_constant_radius))
    
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
    parser = argparse.ArgumentParser()
    parser.add_argument("--length", type=float, required=True, help="square side length")
    parser.add_argument("--points", type=int, required=True, help="number of points in x/y")
    parser.add_argument("--exclude_corners", action='store_true', default=False,
                        help="exclude corner points")
    args = parser.parse_args()
    output_square_discretization_points(args.length,
                                        args.points,
                                        args.exclude_corners)
    
    # if (len(sys.argv) != 5):
    #     print("square_functions.py [length num_points_xy num_neighbors_basis num_neighbors_weight]")
    #     sys.exit()
    # length = float(sys.argv[1])
    # num_points_xy = int(sys.argv[2])
    # num_neighbors_basis = int(sys.argv[3])
    # num_neighbors_weight = int(sys.argv[4])
    
    # output_square_discretization(length,
    #                              num_points_xy,
    #                              num_neighbors_basis,
    #                              num_neighbors_weight,
    #                              True)
    

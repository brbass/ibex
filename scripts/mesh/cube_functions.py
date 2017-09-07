from mesh_functions import *

# Get points and find connectivity for pincell
def get_cube_discretization(length,
                            num_points_xyz,
                            num_neighbors_basis,
                            num_neighbors_weight,
                            use_constant_radius):

    num_points, points = get_cartesian_points_3d(length,
                                                 num_points_xyz)
    
    return get_connectivity(num_points,
                            points,
                            num_neighbors_basis,
                            num_neighbors_weight,
                            use_constant_radius)

def output_cube_discretization_points(length,
                                      num_points_xyz):
    output_path = "cube_{}_{}.xml".format(length, num_points_xyz)
    num_points, points = get_cartesian_points_3d(length,
                                                 num_points_xyz)
    node = xml_points(3, # dimension
                      num_points,
                      points)
    spatial_node = node.find("spatial_discretization")
    
    et.ElementTree(node).write(output_path,
                               pretty_print=True,
                               xml_declaration=True)
    
if __name__ == '__main__':
    if (len(sys.argv) != 3):
        print("cube_functions.py [length num_points_xyz]")
        sys.exit()
    length = float(sys.argv[1])
    num_points_xyz = int(sys.argv[2])
    
    output_cube_discretization_points(length,
                                      num_points_xyz)
    

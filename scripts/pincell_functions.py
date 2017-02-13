from mesh_functions import *

# Get points and find connectivity for pincell
def get_pincell_discretization(radius,
                               length,
                               num_points_xy,
                               num_points_r,
                               num_neighbors_basis,
                               num_neighbors_weight):

    num_points, points = get_pincell_points(radius,
                                            length,
                                            num_points_xy,
                                            num_points_r)
    
    return get_connectivity(num_points,
                            points,
                            num_neighbors_basis,
                            num_neighbors_weight)

def output_pincell_discretization(radius,
                                  length,
                                  num_points_xy,
                                  num_points_r,
                                  num_neighbors_basis,
                                  num_neighbors_weight):
    num_points, points, radius_basis, radius_weight, b_neighbors, b_neighbor_distances, w_neighbors, w_neighbor_distances, bw_neighbors, bw_neighbor_distances \
        = get_pincell_discretization(radius,
                                     length,
                                     num_points_xy,
                                     num_points_r,
                                     num_neighbors_basis,
                                     num_neighbors_weight)
    
    output_path = "in/pincell_{}_{}_{}_{}_{}_{}.xml".format(radius, length, num_points_xy, num_points_r, num_neighbors_basis, num_neighbors_weight)

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
    origins = np.array([[-l2, 0],
                        [l2, 0],
                        [0, -l2],
                        [0, l2]])
    normals = np.array([[-1, 0],
                        [1, 0],
                        [0, -1],
                        [0, 1]])
    boundary_sources = [0, 0, 0, 0]
    for i in range(4):
        surface = et.SubElement(surfaces, "surface")
        surface.set("index", str(i))
        et.SubElement(surface, "shape").text = "plane"
        et.SubElement(surface, "type").text = "boundary"
        et.SubElement(surface, "origin").text = numpy_to_text(origins[i])
        et.SubElement(surface, "normal").text = numpy_to_text(normals[i])
        et.SubElement(surface, "boundary_source").text = str(boundary_sources[i])
    surface = et.SubElement(surfaces, "surface")
    surface.set("index", "4")
    et.SubElement(surface, "shape").text = "cylinder"
    et.SubElement(surface, "type").text = "internal"
    et.SubElement(surface, "origin").text = "0 0"
    et.SubElement(surface, "radius").text = numpy_to_text(np.array(radius))

    regions = et.SubElement(solid, "regions")
    
    region = et.SubElement(regions, "region")
    region.set("index", "0")
    et.SubElement(region, "material").text = str(0)
    sr = et.SubElement(region, "surface_relation")
    sr.set("index", "0")
    et.SubElement(sr, "surface").text = str(4)
    et.SubElement(sr, "relation").text = "inside"

    region = et.SubElement(regions, "region")
    region.set("index", "1")
    et.SubElement(region, "material").text = str(1)
    for i in range(5):
        sr = et.SubElement(region, "surface_relation")
        sr.set("index", str(i))
        et.SubElement(sr, "surface").text = str(i)
        if i < 4:
            et.SubElement(sr, "relation").text = "inside"
        else:
            et.SubElement(sr, "relation").text = "outside"

    et.ElementTree(node).write(output_path,
                               pretty_print=True,
                               xml_declaration=True)

    
def plot_pincell_discretization(radius,
                                length,
                                num_points_xy,
                                num_points_r,
                                num_neighbors_basis,
                                num_neighbors_weight):
    num_points, points, radius_basis, radius_weight, b_neighbors, b_neighbor_distances, w_neighbors, w_neigbor_distances, bw_neighbors, bw_neighbor_distances \
        = get_pincell_discretization(radius,
                                     length,
                                     num_points_xy,
                                     num_points_r,
                                     num_neighbors_basis,
                                     num_neighbors_weight)
    
    k = .7
    range_plot = np.linspace(-length/2, length/2, 50)
    bound_min = -length/2 * np.ones(50)
    bound_max = length/2 * np.ones(50)
    theta = np.linspace(0, 2*np.pi, 50)
    x_rad_plot = radius * np.cos(theta)
    y_rad_plot = radius * np.sin(theta)
    bound_color='k'
    point_color='r'
    support_color='b'
    bound_width = 2.
    for include_radii in [True, False]:
        for d_index, distance in enumerate([radius_basis, radius_weight]):
            if d_index == 0:
                dist_type = "basis"
                num_neighbors = num_neighbors_basis
            else:
                dist_type = "weight"
                num_neighbors = num_neighbors_weight
            plt.figure()
            plt.axis('equal')
            plt.xlim(-k*length, k*length)
            plt.ylim(-k*length, k*length)
            plt.xlabel("x")
            plt.ylabel("y")
            if include_radii:
                for i in range(num_points):
                    rmax = 1*distance[i]
                    x_vals = points[i, 0] + rmax * np.cos(theta)
                    y_vals = points[i, 1] + rmax * np.sin(theta)
                    plt.plot(x_vals, y_vals, color=support_color)
                plt.plot(x_vals, y_vals, color=support_color, label="support")
            plt.plot(range_plot, bound_min, color=bound_color, linewidth=bound_width, label="boundaries")
            plt.plot(range_plot, bound_max, color=bound_color, linewidth=bound_width)
            plt.plot(bound_min, range_plot, color=bound_color, linewidth=bound_width)
            plt.plot(bound_max, range_plot, color=bound_color, linewidth=bound_width)
            plt.plot(x_rad_plot, y_rad_plot, color=bound_color, linewidth=bound_width)
            plt.scatter(points[:,0], points[:,1], color=point_color, marker='o', label="centers")
            plt.legend()
            plt.savefig("figs/rbf_disc_{}_{}_{}_{}_{}_{}_{}.pdf".format(radius, length, num_points_xy, num_points_r, num_neighbors, dist_type, include_radii), bbox_inches='tight')
            plt.close()

if __name__ == '__main__':
    if (len(sys.argv) != 7):
        print("pincell_functions.py [radius length num_points_xy num_points_r num_neighbors_basis num_neighbors_weight]")
        sys.exit()
    radius = float(sys.argv[1])
    length = float(sys.argv[2])
    num_points_xy = int(sys.argv[3])
    num_points_r = int(sys.argv[4])
    num_neighbors_basis = int(sys.argv[5])
    num_neighbors_weight = int(sys.argv[6])
    
    output_pincell_discretization(radius,
                                  length,
                                  num_points_xy,
                                  num_points_r,
                                  num_neighbors_basis,
                                  num_neighbors_weight)
    

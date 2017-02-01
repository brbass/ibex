import sys, io
import numpy as np
import scipy as sp
import scipy.spatial as sps
from matplotlib import pyplot as plt
import lxml.etree as et

# Get a Cartesian mesh with a cylindrical cutout
def get_cartesian_mesh_cutout(radius, # Radius of cutout circle
                              length, # Length of a single side
                              num_points_xy, # Number of points in one dimension
                              delta): # Extra distance outside of radius to avoid putting points
    # Get total number of points and the actual points
    num_points = num_points_xy**2
    points_x = np.linspace(-length/2, length/2, num_points_xy, endpoint=True)
    points_y = np.linspace(-length/2, length/2, num_points_xy, endpoint=True)

    # Create array of points outside cutout circle
    points = np.zeros((num_points, 2))
    point = 0
    for i, x in enumerate(points_x):
        for j, y in enumerate(points_y):
            if x**2 + y**2 > (delta + radius)**2:
                points[point, :] = [x, y]
                point += 1
    points = points[0:point]
    return points

# Get a cylindrical mesh
def get_cylindrical_points(radius, # Radius of mesh
                           num_points_r, # Number of points in radial direction
                           radius_endpoint = False): # Include the outer edge of the circle
    # Get radius values
    radius_vals = np.linspace(0, radius, num_points_r, endpoint=radius_endpoint)
    dr = radius_vals[1] - radius_vals[0]
    max_points = 2*int(np.power(2 * radius / dr, 2))

    # Get points
    points = np.zeros((max_points, 2))
    point = 0
    for r in radius_vals:
        if r == 0: # Center point
            points[point, :] = [0, 0]
            point += 1
        else: # Place points in a circle of constant radius with distance approximately equal to dr
            # Get angular values
            dt = np.arcsin(dr / r)
            num_theta = int(np.ceil(2 * np.pi / dt))
            if point == 1:
                num_theta *= 2
            theta_vals = np.linspace(0, 2 * np.pi, num_theta, endpoint=False)
            # Convert values to x-y coordinates
            for t in theta_vals:
                points[point, :] = [r * np.cos(t), r * np.sin(t)]
                point += 1
    points = points[0:point]

    return points

# Combine cylindrical and Cartesian meshes
def get_pincell_points(radius,
                       length,
                       num_points_xy,
                       num_points_r):
    cyl_points = get_cylindrical_points(radius,
                                        num_points_r,
                                        True)
    cart_points = get_cartesian_mesh_cutout(radius,
                                            length,
                                            num_points_xy,
                                            0.7 * radius / num_points_r)
    points = np.concatenate((cart_points, cyl_points), axis=0)
    num_points = len(points)

    return num_points, points

# Get kd_tree and neighbor information
def get_neighbor_information(num_points, # Total number of points
                             num_neighbors, # Desired number of points contained within the calculated radius
                             points): # Locations of points
    # Get kd-tree
    kd_tree = sps.KDTree(points)

    support_radius = np.zeros(num_points)
    neighbors = np.zeros((num_points, num_neighbors))
    neighbor_distances = np.zeros((num_points, num_neighbors))
    for i, point in enumerate(points):
        # Get distance that contains num_neighbors other points within radius and the indices of the other points
        d, indices = kd_tree.query(point, k=num_neighbors)
        neighbor_distances[i, :] = d
        neighbors[i, :] = indices
        
    return neighbors, neighbor_distances, kd_tree

# Get neighbors
def get_neighbors(num_points,
                  points,
                  neighbor_pairs, # List of possible pairs of neighbors
                  radius_base, # Max distance for the desired basis
                  radius_other): # Max distance for the other set to compare against
    neighbors = []
    neighbor_distances = []
    num_neighbors = np.zeros(num_points, dtype=int)
    for i, local_neighbors in enumerate(neighbor_pairs):
        base_neighbors = []
        base_distances = []
        for j in local_neighbors:
            dist = np.sqrt(np.sum(np.power(points[i] - points[j], 2)))
            rad = radius_other[j] + radius_base[i]
            if dist < rad:
                base_neighbors.append(j)
                base_distances.append(dist)
        base_neighbors = np.asarray(base_neighbors)
        base_distances = np.asarray(base_distances)
        indices = np.argsort(base_distances)
        base_neighbors = base_neighbors[indices]
        base_distances = base_distances[indices]
        neighbors.append(base_neighbors)
        neighbor_distances.append(base_distances)
        num_neighbors[i] = len(base_neighbors)
        
    neighbors = np.asarray(neighbors)
    neighbor_distances = np.asarray(neighbor_distances)

    return num_neighbors, neighbors, neighbor_distances
    
# Get connectivity for points
def get_connectivity(num_points,
                     points,
                     num_neighbors_basis,
                     num_neighbors_weight):
    # Get information on point connectivity
    num_neighbors = max(num_neighbors_basis, num_neighbors_weight) + 1
    temp_neighbors, temp_neighbor_distances, kd_tree \
        = get_neighbor_information(num_points,
                                   num_neighbors,
                                   points)

    # Get the max distance for each of the basis and weight functions
    mult = 1.
    radius_basis = mult * temp_neighbor_distances[:, num_neighbors_basis]
    radius_weight = mult * temp_neighbor_distances[:, num_neighbors_weight]
    
    # Get the maximum radius for each of the basis and weight functions
    overall_max_basis = np.amax(radius_basis)
    overall_max_weight = np.amax(radius_weight)
    search_radius = overall_max_basis + overall_max_weight

    # Get pairs of points that are less than the maximum distance apart
    neighbor_pairs = kd_tree.query_ball_tree(kd_tree, search_radius)

    # Get basis, weight and basis/weight functions that intersect
    b_num_neighbors, b_neighbors, b_neighbor_distances \
        = get_neighbors(num_points,
                        points,
                        neighbor_pairs,
                        radius_basis,
                        radius_basis)
    w_num_neighbors, w_neighbors, w_neighbor_distances \
        = get_neighbors(num_points,
                        points,
                        neighbor_pairs,
                        radius_weight,
                        radius_weight)
    bw_num_neighbors, bw_neighbors, bw_neighbor_distances \
        = get_neighbors(num_points,
                        points,
                        neighbor_pairs,
                        radius_weight,
                        radius_basis)
    
    return num_points, points, radius_basis, radius_weight, b_neighbors, b_neighbor_distances, w_neighbors, w_neighbor_distances, bw_neighbors, bw_neighbor_distances

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

# Convert data to text
def numpy_to_text(data):
    a = np.array2string(data, separator=" ", precision=20, max_line_width=1e5)
    if a[0] == "[":
        return a[1:-1]
    else:
        return a

# Get XML node describing discretization
def xml_discretization(output_path,
                       num_points,
                       points,
                       radius_basis,
                       radius_weight,
                       b_neighbors,
                       b_neighbor_distances,
                       w_neighbors,
                       w_neighbor_distances,
                       bw_neighbors,
                       bw_neighbor_distances):
    node = et.Element("discretization")
    et.SubElement(node, "number_of_points").text = str(num_points)
    
    basis_top = et.SubElement(node, "basis_functions")
    for i in range(num_points):
        basis = et.SubElement(basis_top, "basis")
        basis.set("index", str(i))
        num_b_neighbors = len(b_neighbors[i])
        et.SubElement(basis, "radius").text = numpy_to_text(radius_basis[i])
        et.SubElement(basis, "position").text = numpy_to_text(points[i, :])
        et.SubElement(basis, "number_of_basis_neighbors").text = str(num_b_neighbors)
        et.SubElement(basis, "basis_neighbors").text = numpy_to_text(b_neighbors[i])
        et.SubElement(basis, "basis_neighbor_distances").text = numpy_to_text(b_neighbor_distances[i])
        
    weight_top = et.SubElement(node, "weight_functions")
    for i in range(num_points):
        weight = et.SubElement(weight_top, "weight")
        weight.set("index", str(i))
        num_b_neighbors = len(b_neighbors[i])
        num_bw_neighbors = len(bw_neighbors[i])
        et.SubElement(weight, "radius").text = numpy_to_text(radius_weight[i])
        et.SubElement(weight, "position").text = numpy_to_text(points[i, :])
        et.SubElement(weight, "number_of_basis_neighbors").text = str(num_b_neighbors)
        et.SubElement(weight, "basis_neighbors").text = numpy_to_text(bw_neighbors[i])
        et.SubElement(weight, "basis_neighbor_distances").text = numpy_to_text(bw_neighbor_distances[i])
        et.SubElement(weight, "number_of_weight_neighbors").text = str(num_bw_neighbors)
        et.SubElement(weight, "weight_neighbors").text = numpy_to_text(w_neighbors[i])
        et.SubElement(weight, "weight_neighbor_distances").text = numpy_to_text(w_neighbor_distances[i])
    return node

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
    
    output_path = "in/rbf_disc_{}_{}_{}_{}_{}_{}.xml".format(radius, length, num_points_xy, num_points_r, num_neighbors_basis, num_neighbors_weight)

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
    
    et.ElementTree(node).write(output_path, pretty_print=True, xml_declaration=True)

    
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
        print("mesh_functions.py [radius length num_points_xy num_points_r num_neighbors_basis num_neighbors_weight]")
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
        
    
                       


import sys, io
import numpy as np
import scipy as sp
import scipy.spatial as sps
from matplotlib import pyplot as plt
import lxml.etree as et

# Get a Cartesian mesh with a cylindrical cutout
def get_cartesian_mesh_cutout(radius,
                              length,
                              num_points_xy,
                              delta):
    num_points = num_points_xy**2
    points_x = np.linspace(-length/2, length/2, num_points_xy, endpoint=True)
    points_y = np.linspace(-length/2, length/2, num_points_xy, endpoint=True)
    
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
def get_cylindrical_points(radius,
                           num_points_r,
                           radius_endpoint = False):
    radius_vals = np.linspace(0, radius, num_points_r, endpoint=radius_endpoint)
    dr = radius_vals[1] - radius_vals[0]
    max_points = 2*int(np.power(2 * radius / dr, 2))
    points = np.zeros((max_points, 2))
    point = 0
    for r in radius_vals:
        if r == 0:
            points[point, :] = [0, 0]
            point += 1
        else:
            dt = np.arcsin(dr / r)
            num_theta = int(np.ceil(2 * np.pi / dt))
            if point == 1:
                num_theta *= 2
            theta_vals = np.linspace(0, 2 * np.pi, num_theta, endpoint=False)
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
def get_neighbor_information(num_points,
                             num_neighbors,
                             points):
    kd_tree = sps.KDTree(points)
    
    support_radius = np.zeros(num_points)
    neighbors = np.zeros((num_points, num_neighbors))
    neighbor_distances = np.zeros((num_points, num_neighbors))
    for i, point in enumerate(points):
        d, indices = kd_tree.query(point, k=num_neighbors)
        neighbor_distances[i, :] = d
        neighbors[i, :] = indices

    return neighbors, neighbor_distances, kd_tree

# Get connectivity for points
def get_connectivity(num_points,
                     points,
                     num_neighbors_basis,
                     num_neighbors_weight):
    num_neighbors = max(num_neighbors_basis, num_neighbors_weight) + 1
    temp_neighbors, temp_neighbor_distances, kd_tree \
        = get_neighbor_information(num_points,
                                   num_neighbors,
                                   points)
    
    mult = 1.
    max_distance_basis = mult * temp_neighbor_distances[:, num_neighbors_basis]
    max_distance_weight = mult * temp_neighbor_distances[:, num_neighbors_weight]
    
    overall_max_basis = np.amax(max_distance_basis)
    overall_max_weight = np.amax(max_distance_weight)
    search_radius = overall_max_basis + overall_max_weight
    neighbor_pairs = kd_tree.query_ball_tree(kd_tree, search_radius)
    
    neighbors = []
    neighbor_distances = []
    num_neighbors = np.zeros(num_points, dtype=int)
    for i, local_neighbors in enumerate(neighbor_pairs):
        weight_neighbors = []
        weight_distances = []
        for j in local_neighbors:
            dist = np.sqrt(np.sum(np.power(points[i] - points[j], 2)))
            rad = max_distance_basis[j] + max_distance_weight[i]
            if dist < rad:
                weight_neighbors.append(j)
                weight_distances.append(dist)
        weight_neighbors = np.asarray(weight_neighbors)
        weight_distances = np.asarray(weight_distances)
        indices = np.argsort(weight_distances)
        weight_neighbors = weight_neighbors[indices]
        weight_distances = weight_distances[indices]
        neighbors.append(weight_neighbors)
        neighbor_distances.append(weight_distances)
        num_neighbors[i] = len(weight_neighbors)
        
    neighbors = np.asarray(neighbors)
    neighbor_distances = np.asarray(neighbor_distances)
    
    return num_points, points, max_distance_basis, max_distance_weight, neighbors, neighbor_distances

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

def numpy_to_text(data):
    a = np.array2string(data, separator=" ", precision=20, max_line_width=1e5)
    if a[0] == "[":
        return a[1:-1]
    else:
        return a

def xml_discretization(output_path,
                       num_points,
                       points,
                       max_distance_basis,
                       max_distance_weight,
                       neighbors,
                       neighbor_distances):
    node = et.Element("discretization")
    et.SubElement(node, "number_of_points").text = str(num_points)
    
    basis_top = et.SubElement(node, "basis_functions")
    for i in range(num_points):
        basis = et.SubElement(basis_top, "basis")
        basis.set("index", str(i))
        et.SubElement(basis, "max_distance").text = numpy_to_text(max_distance_basis[i])
        et.SubElement(basis, "position").text = numpy_to_text(points[i, :])
        
    weight_top = et.SubElement(node, "weight_functions")
    for i in range(num_points):
        weight = et.SubElement(weight_top, "weight")
        weight.set("index", str(i))
        num_neighbors = len(neighbors[i])
        et.SubElement(weight, "number_of_neighbors").text = str(num_neighbors)
        et.SubElement(weight, "max_distance").text = numpy_to_text(max_distance_weight[i])
        et.SubElement(weight, "position").text = numpy_to_text(points[i, :])
        et.SubElement(weight, "neighbors").text = numpy_to_text(neighbors[i])
        et.SubElement(weight, "neighbor_distances").text = numpy_to_text(neighbor_distances[i])
    return node

def output_pincell_discretization(radius,
                                  length,
                                  num_points_xy,
                                  num_points_r,
                                  num_neighbors_basis,
                                  num_neighbors_weight):
    num_points, points, max_distance_basis, max_distance_weight, neighbors, neighbor_distances \
        = get_pincell_discretization(radius,
                                     length,
                                     num_points_xy,
                                     num_points_r,
                                     num_neighbors_basis,
                                     num_neighbors_weight)
    
    output_path = "in/rbf_disc_{}_{}_{}_{}_{}_{}.xml".format(radius, length, num_points_xy, num_points_r, num_neighbors_basis, num_neighbors_weight)
    
    node = xml_discretization(output_path,
                              num_points,
                              points,
                              max_distance_basis,
                              max_distance_weight,
                              neighbors,
                              neighbor_distances)

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
    et.SubElement(region, "material").text = str(0)
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
    num_points, points, max_distance_basis, max_distance_weight, neighbors, neighbor_distances \
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
        for d_index, distance in enumerate([max_distance_basis, max_distance_weight]):
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
        
    
                       


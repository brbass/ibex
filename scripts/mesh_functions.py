import sys, io
import numpy as np
import scipy as sp
import scipy.spatial as sps
from matplotlib import pyplot as plt
import lxml.etree as et

# Get a 1D Cartesian mesh
def get_cartesian_points_1d(xmin,
                            xmax,
                            num_points):
    points_1d = np.linspace(xmin, xmax, num_points, endpoint=True)
    points = np.zeros((num_points, 1))
    for i in range(num_points):
        points[i][0] = points_1d[i]
    return points

# Get a 2D Cartesian mesh
def get_cartesian_points_2d(xmin,
                            xmax,
                            ymin,
                            ymax,
                            numx,
                            numy,
                            exclude_corners = False):
    # Get total number of points and the actual points
    num_points = numx * numy;
    points_x = np.linspace(xmin, xmax, numx, endpoint=True)
    points_y = np.linspace(ymin, ymax, numy, endpoint=True)
    
    # Create array of points outside cutout circle
    points = np.zeros((num_points, 2))
    point = 0
    for i, x in enumerate(points_x):
        for j, y in enumerate(points_y):
                points[i + numx * j, :] = [x, y]
    if exclude_corners:
        xc = numx - 1
        yc = numy - 1
        corner_indices = [0, xc, numx * yc, xc + numx * yc]
        points = np.delete(points, corner_indices, axis=0)
        num_points = num_points - 4
        
    return num_points, points
    
# Get a 3D Cartesian mesh
def get_cartesian_points_3d(xmin,
                            xmax,
                            ymin,
                            ymax,
                            zmin,
                            zmax,
                            numx,
                            numy,
                            numz,
                            exclude_corners = False):
    # Get total number of points and the actual points
    num_points = numx * numy * numz
    points_x = np.linspace(xmin, xmax, numx, endpoint=True)
    points_y = np.linspace(ymin, ymax, numy, endpoint=True)
    points_z = np.linspace(zmin, zmax, numz, endpoint=True)
    
    # Create array of points outside cutout circle
    points = np.zeros((num_points, 3))
    point = 0
    for i, x in enumerate(points_x):
        for j, y in enumerate(points_y):
            for k, z in enumerate(points_z):
                points[i + numx * (j + numy * k), :] = [x, y, z]
    if exclude_corners:
        corner_indices = []
        xc = numx - 1
        yc = numy - 1
        zc = numz - 1
        for i, x in enumerate(points_x):
            corner_indices.append(i + numx * (0 + numy * 0))
            corner_indices.append(i + numx * (yc + numy * 0))
            corner_indices.append(i + numx * (0 + numy * zc))
            corner_indices.append(i + numx * (yc + numy * zc))
        for j, y in enumerate(points_y):
            corner_indices.append(0 + numx * (j + numy * 0)) 
            corner_indices.append(xc + numx * (j + numy * 0))
            corner_indices.append(0 + numx * (j + numy * zc))
            corner_indices.append(xc + numx * (j + numy * zc))
        for k, z in enumerate(points_z):
            corner_indices.append(0 + numx * (0 + numy * k)) 
            corner_indices.append(xc + numx * (0 + numy * k))
            corner_indices.append(0 + numx * (yc + numy * k))
            corner_indices.append(xc + numx * (yc + numy * k))
        corner_indices = np.unique(corner_indices)
        points = np.delete(points, corner_indices, axis=0)
        num_points = num_points - len(corner_indices)
    return num_points, points
    
# Get slab with automatic min/max
def get_slab(length,
             num_points):
    points_1d = np.linspace(-length/2, length/2, num_points, endpoint=True)
    points = np.zeros((num_points, 1))
    for i in range(num_points):
        points[i][0] = points_1d[i]
    return points

# Get a square with automatic min/max
def get_square(length, # Length of a single side
               num_points_xy, # Number of points in one dimension
               exclude_corners = False): # Exclude points in corners
    # Get total number of points and the actual points
    num_points = num_points_xy**2
    points_x = np.linspace(-length/2, length/2, num_points_xy, endpoint=True)
    points_y = np.linspace(-length/2, length/2, num_points_xy, endpoint=True)
    
    # Create array of points outside cutout circle
    points = np.zeros((num_points, 2))
    point = 0
    for i, x in enumerate(points_x):
        for j, y in enumerate(points_y):
                points[i + num_points_xy * j, :] = [x, y]
    if exclude_corners:
        k = num_points_xy - 1
        corner_indices = [0, k, 0 + num_points_xy * k, k + num_points_xy * k]
        points = np.delete(points, corner_indices, axis=0)
        num_points = num_points - 4
    return num_points, points
    
# Get a Cartesian mesh with a cylindrical cutout
def get_square_cutout(radius, # Radius of cutout circle
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

# Get a 3D Cartesian mesh
def get_cube(length, # Length of a single side
             num_points_xyz): # Number of points in one dimension
    # Get total number of points and the actual points
    num_points = num_points_xyz**3
    points_x = np.linspace(-length/2, length/2, num_points_xyz, endpoint=True)
    points_y = np.linspace(-length/2, length/2, num_points_xyz, endpoint=True)
    points_z = np.linspace(-length/2, length/2, num_points_xyz, endpoint=True)
    
    # Create array of points outside cutout circle
    points = np.zeros((num_points, 3))
    point = 0
    for i, x in enumerate(points_x):
        for j, y in enumerate(points_y):
            for k, z in enumerate(points_z):
                points[i + num_points_xyz * (j + num_points_xyz * k), :] = [x, y, z]
    return num_points, points

# Get a cylindrical mesh
def get_cylindrical_points(radius, # Radius of mesh
                           num_points_r, # Number of points in radial direction
                           include_endpoint = False): # Include the outer edge of the circle
    return get_annular_points(0,
                              radius,
                              num_points_r,
                              include_endpoint)

# Get a ring of points
def get_ring_points(r,
                    dr,
                    randomize_start=False):
    if r == 0.0:
        points = [[0, 0]]
    else:
        if dr > r:
            dr = r
        dt = np.arcsin(dr / r)
        num_theta = int(np.ceil(2 * np.pi / dt))
        if num_theta < 4:
            num_theta = 4
        
        if randomize_start:
            theta0 = np.random.rand() * 2 * np.pi
            theta1 = theta0 + 2 * np.pi
        else:
            theta0 = 0.
            theta1 = 2 * np.pi
        theta_vals = np.linspace(theta0, theta1, num_theta, endpoint=False)
        points = [[r * np.cos(t), r * np.sin(t)] for t in theta_vals]

    return points

# Get annular points
def get_annular_points(r1, # inside radius
                       r2, # outside radius
                       num_points_r,
                       include_endpoint = True):
    # Get radius values
    radius_vals = np.linspace(r1, r2, num_points_r, endpoint=include_endpoint)
    if include_endpoint:
        dr = radius_vals[1] - radius_vals[0]
    else:
        dr = (r2 - r1) / num_points_r
    max_points = 2*int(np.power(2 * r2 / dr, 2))
    
    # Get points
    points = np.zeros((max_points, 2))
    point = 0
    for r in radius_vals:
        points_temp = get_ring_points(r,
                                      dr)
        for point_temp in points_temp:
            points[point, :] = point_temp
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
    cart_points = get_square_cutout(radius,
                                    length,
                                    num_points_xy,
                                    0.7 * radius / (num_points_r - 1))
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
            rad_delta = rad * (1. - 1e-10)
            if dist < rad_delta:
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
                     num_neighbors_weight,
                     use_constant_radius = False):
    # Get information on point connectivity
    num_neighbors = max(num_neighbors_basis, num_neighbors_weight) + 1
    temp_neighbors, temp_neighbor_distances, kd_tree \
        = get_neighbor_information(num_points,
                                   num_neighbors,
                                   points)

    # Get the max distance for each of the basis and weight functions
    mult = 1.
    if use_constant_radius:
        radius_basis = mult * temp_neighbor_distances[:, 1] * 5
        radius_weight = mult * temp_neighbor_distances[:, 1] * 5
    else:
        radius_basis = mult * temp_neighbor_distances[:, num_neighbors_basis]
        radius_weight = mult * temp_neighbor_distances[:, num_neighbors_weight]
    
    # Get the maximum radius for each of the basis and weight functions
    overall_max_basis = np.amax(radius_basis)
    overall_max_weight = np.amax(radius_weight)
    overall_min_basis = np.amax(radius_basis)
    overall_min_weight = np.amin(radius_weight)
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

# Convert data to text
def numpy_to_text(data):
    np.set_printoptions(threshold=np.inf)
    a = np.array2string(data, separator=" ", precision=16, max_line_width=1e8)
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
    node = et.Element("spatial_discretization")
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
        # et.SubElement(basis, "basis_neighbor_distances").text = numpy_to_text(b_neighbor_distances[i])
        
    weight_top = et.SubElement(node, "weight_functions")
    for i in range(num_points):
        weight = et.SubElement(weight_top, "weight")
        weight.set("index", str(i))
        num_w_neighbors = len(w_neighbors[i])
        num_bw_neighbors = len(bw_neighbors[i])
        et.SubElement(weight, "radius").text = numpy_to_text(radius_weight[i])
        et.SubElement(weight, "position").text = numpy_to_text(points[i, :])
        et.SubElement(weight, "number_of_basis_neighbors").text = str(num_bw_neighbors)
        et.SubElement(weight, "basis_neighbors").text = numpy_to_text(bw_neighbors[i])
        # et.SubElement(weight, "basis_neighbor_distances").text = numpy_to_text(bw_neighbor_distances[i])
        et.SubElement(weight, "number_of_weight_neighbors").text = str(num_w_neighbors)
        et.SubElement(weight, "weight_neighbors").text = numpy_to_text(w_neighbors[i])
        # et.SubElement(weight, "weight_neighbor_distances").text = numpy_to_text(w_neighbor_distances[i])
    return node

def xml_points(dimension,
               num_points,
               points):
    flat_points = points.flatten()
    
    node = et.Element("input")
    spatial_node = et.SubElement(node, "spatial_discretization")
    et.SubElement(spatial_node, "dimension").text = str(dimension)
    et.SubElement(spatial_node, "number_of_points").text = str(num_points)
    et.SubElement(spatial_node, "points").text = numpy_to_text(flat_points)
    
    return node
    

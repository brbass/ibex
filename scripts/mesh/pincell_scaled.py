from mesh_functions import *
import numpy as np
import argparse

def get_points(num_radii,
               max_delta,
               initial_delta,
               ring_start = 0.4101,
               length = 1.26,
               randomize_start=False,
               exclude_corners = False):
    # Set parameters
    initial_delta = np.abs(initial_delta)
    r0 = ring_start - initial_delta / 2
    r1 = ring_start + initial_delta / 2
    r = 0.0
    dr0 = -initial_delta
    a = np.power((r - r0) / dr0, 1. / (num_radii - 1.))
    l = length
    l2 = l / 2
    print("delta multiplication factor: a = ", a)
    ind_add = int(np.ceil(l2 / max_delta))
    
    # Get points and spacing for each radius
    points_ifba = [r0, r1]
    delta_ifba = [initial_delta, initial_delta]
    delta_fuel = [min(max_delta, initial_delta * np.power(a, n + 1)) for n in range(0, num_radii + ind_add)]
    points_fuel = [r0 - delta_fuel[0]]
    for i in range(num_radii + ind_add - 1):
        point = points_fuel[-1] - delta_fuel[i + 1]
        if point >= 0:
            if point < delta_fuel[-2] / 4:
                point = 0.
            points_fuel.append(point)
    delta_fuel = delta_fuel[0:len(points_fuel)]
    delta_mod = [min(max_delta, initial_delta * np.power(a, n + 1)) for n in range(0, num_radii + ind_add)]
    points_mod = [r1 + delta_mod[0]]
    for i in range(num_radii + ind_add - 1):
        point = points_mod[-1] + delta_mod[i + 1]
        points_mod.append(point)
    
    # Combine sets of points
    points_r = np.concatenate((points_fuel, points_ifba, points_mod))
    delta_r = np.abs(np.concatenate((delta_fuel, delta_ifba, delta_mod)))
    indices = np.argsort(points_r)
    points_r = points_r[indices]
    delta_r = delta_r[indices]
    spacing_r = [delta_r[0]]
    for i in range(1, len(points_r) - 1):
        spacing_r.append((delta_r[i + 1] + delta_r[i - 1]) / 2)
    spacing_r.append(delta_r[-1])

    # Get full set of points
    points = []
    for i, (r, spacing) in enumerate(zip(points_r, spacing_r)):
        points_temp = get_ring_points(r, spacing, randomize_start)
        for point_temp in points_temp:
            points.append(point_temp)
            
    # Get distance to ignore near boundary
    index = (np.abs(points_r - l2)).argmin()
    delta_l = delta_r[index]
    l_ignore = (l2 - 0.7 * delta_l)

    # Remove those points outside boundary
    points_all = points
    points = []
    for point in points_all:
        if np.abs(point[0]) < l_ignore and np.abs(point[1]) < l_ignore:
            points.append(point)

    # Add points along boundary
    num_intervals_l = int(np.ceil(l / delta_l))
    num_points_l = num_intervals_l + 1
    delta_l = l / num_intervals_l
    for sign in [-1, 1]:
        if exclude_corners:
            for i in range(1, num_points_l - 1):
                pos = -l2 + delta_l * i
                points.append([pos, sign * l2])
        else:
            for i in range(num_points_l):
                pos = -l2 + delta_l * i
                points.append([pos, sign * l2])
        for i in range(1, num_points_l - 1):
            pos = -l2 + delta_l * i
            points.append([sign * l2, pos])
            
    # Return result
    num_points = len(points)
    points = np.array(points)
    return num_points, points

def output_points(num_radii,
                  max_delta,
                  initial_delta,
                  ring_start = 0.4101,
                  length = 1.26,
                  randomize_start = True,
                  exclude_corners = False,
                  restriction = None,
                  restrict_circle = None,
                  plot = False):
    # Get points and output path
    num_points, points = get_points(num_radii, max_delta, initial_delta,
                                    length=length,
                                    ring_start=ring_start,
                                    randomize_start=randomize_start,
                                    exclude_corners=exclude_corners)
    if restriction != None:
        old_points = points
        points = []
        x1 = restriction[0]
        x2 = restriction[1]
        y1 = restriction[2]
        y2 = restriction[3]
        for point in old_points:
            if x1 <= point[0] and point[0] <= x2 and y1 <= point[1] and point[1] <= y2:
                points.append(point)
        num_points_inside = len(points)
        num_points_l = int(2.0 * np.sqrt(num_points_inside))
        delta_l_x = (x2 - x1) / (num_points_l - 1)
        delta_l_y = (y2 - y1) / (num_points_l - 1)
        for xval, yval in zip([x1, x2], [y1, y2]):
            for i in range(num_points_l):
                pos = x1 + delta_l_x * i
                points.append([pos, yval])
            for i in range(1, num_points_l - 1):
                pos = y1 + delta_l_y * i
                points.append([xval, pos])
        points = np.array(points)
        num_points = len(points)
    if restrict_circle != None:
        radius = restrict_circle[0]
        delta_frac = restrict_circle[1]
        old_points = points
        points = []
        for point in old_points:
            if point[0] * point[0] + point[1] * point[1] <= radius * radius:
                points.append(point)
        num_points_inside = len(points)
        num_points_l = int(np.sqrt(num_points_inside))
        delta_l = delta_frac * initial_delta
        old_points = points
        points = []
        for point in old_points:
            rmax = radius - 0.7 * delta_l
            if point[0] * point[0] + point[1] * point[1] <= rmax * rmax:
                points.append(point)
        points_temp = get_ring_points(radius, max_delta, randomize_start)
        for point_temp in points_temp:
            points.append(point_temp)
        points = np.array(points)
        num_points = len(points)
    print("mesh has {} points".format(num_points))
    if exclude_corners:
        corner_string = "_exclude"
    else:
        corner_string = ""
    if randomize_start:
        corner_string += "_rand"
    if restrict_circle:
        corner_string += "_circ{}".format(restrict_circle[0])
    output_path = "pincell_scaled_{}_{}_{}_{}_{}{}.xml".format(num_radii,
                                                               max_delta,
                                                               initial_delta,
                                                               ring_start,
                                                               length,
                                                               corner_string)

    # Get node and add input data
    node = xml_points(2, # dimension
                      num_points,
                      points)
    et.SubElement(node, "num_radii").text = str(num_radii)
    et.SubElement(node, "max_delta").text = str(max_delta)
    et.SubElement(node, "initial_delta").text = str(initial_delta)
    
    # Output xml node
    et.ElementTree(node).write(output_path,
                               pretty_print=True,
                               xml_declaration=True)
    
    # Plot if desired
    if plot:
        plt.figure()
        plt.scatter(points[:, 0], points[:, 1], s=2)
        plt.axes().set_aspect('equal')
        plt.show()
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--num_radii", type=int, required=True, help="number of points inside the radius")
    parser.add_argument("--initial_delta", type=float, required=True)
    parser.add_argument("--max_delta", type=float, required=True)
    parser.add_argument("--exclude_corners", action='store_true', default=False,
                        help="exclude corner points")
    parser.add_argument("--randomize_start", action='store_true', default=False,
                        help="randomize ring starting points")
    parser.add_argument("--length", type=float, required=False, default=1.26,
                        help="square side length")
    parser.add_argument("--ring_start", type=float, required=False, default=0.4101,
                        help="center location of the ring")
    parser.add_argument("--restrict", type=float, nargs=4, required=False, default=None,
                        help="restrict points to those from x1,x2,y1,y2")
    parser.add_argument("--restrict_circle", type=float, nargs=2, required=False, default=None,
                        help="restrict points to those within the given radius, args are radius and fraction of initial_delta from edge")
    parser.add_argument("--plot", action='store_true', default=False,
                        help="plot points")
    args = parser.parse_args()
    output_points(args.num_radii,
                  args.max_delta,
                  args.initial_delta,
                  ring_start=args.ring_start,
                  length=args.length,
                  randomize_start=args.randomize_start,
                  exclude_corners=args.exclude_corners,
                  restriction=args.restrict,
                  restrict_circle=args.restrict_circle,
                  plot=args.plot)

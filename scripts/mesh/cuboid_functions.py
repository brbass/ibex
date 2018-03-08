from mesh_functions import *
import argparse

def output_cuboid_points(xmin,
                         xmax,
                         ymin,
                         ymax,
                         zmin,
                         zmax,
                         numx,
                         numy,
                         numz,
                         exclude_corners = False):
    output_path = "cuboid_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}.xml".format(xmin,
                                                                    xmax,
                                                                    ymin,
                                                                    ymax,
                                                                    zmin,
                                                                    zmax,
                                                                    numx,
                                                                    numy,
                                                                    numz,
                                                                    exclude_corners)
    num_points, points = get_cartesian_points_3d(xmin,
                                                 xmax,
                                                 ymin,
                                                 ymax,
                                                 zmin,
                                                 zmax,
                                                 numx,
                                                 numy,
                                                 numz,
                                                 exclude_corners)
    
    node = xml_points(3, # dimension
                      num_points,
                      points)
    spatial_node = node.find("spatial_discretization")
    
    et.ElementTree(node).write(output_path,
                               pretty_print=True,
                               xml_declaration=True)
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--x", type=float, nargs=2, required=True,
                        help="xmin and xmax")
    parser.add_argument("--y", type=float, nargs=2, required=True,
                        help="ymin and ymax")
    parser.add_argument("--z", type=float, nargs=2, required=True,
                        help="zmin and zmax")
    parser.add_argument("--num", type=int, nargs=3, required=True,
                        help="num dimensional points for x, y, z")
    parser.add_argument("--exclude_corners", action='store_true', default=False,
                        help="exclude corner points")
    args = parser.parse_args()
    output_cuboid_points(args.x[0],
                         args.x[1],
                         args.y[0],
                         args.y[1],
                         args.z[0],
                         args.z[1],
                         args.num[0],
                         args.num[1],
                         args.num[2],
                         args.exclude_corners)
    

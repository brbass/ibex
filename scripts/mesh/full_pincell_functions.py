from mesh_functions import *

def get_points_full_pincell(r_fuel,
                            r_ifba,
                            r_gap,
                            r_clad,
                            l_mod,
                            num_points_fuel,
                            num_points_ifba,
                            num_points_gap,
                            num_points_clad,
                            num_points_mod):
    mod_points = get_cartesian_mesh_cutout(r_clad,
                                           l_mod,
                                           num_points_mod,
                                           0.7*(r_clad - r_gap) / num_points_clad)
    fuel_points = get_cylindrical_points(r_fuel,
                                         num_points_fuel,
                                         False)
    ifba_points = get_annular_points(r_fuel,
                                     r_ifba,
                                     num_points_ifba,
                                     False)
    gap_points = get_annular_points(r_ifba,
                                    r_gap,
                                    num_points_gap,
                                    False)
    clad_points = get_annular_points(r_gap,
                                     r_clad,
                                     num_points_clad,
                                     True)
    
    points = np.concatenate((mod_points, clad_points, gap_points, ifba_points, fuel_points), axis=0)
    num_points = len(points)
    
    return num_points, points

def output_full_pincell(r_fuel,
                        r_ifba,
                        r_gap,
                        r_clad,
                        l_mod,
                        num_points_fuel,
                        num_points_ifba,
                        num_points_gap,
                        num_points_clad,
                        num_points_mod):
    output_path = "in/full_pincell_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}.xml".format(r_fuel,
                                                                             r_ifba,
                                                                             r_gap,
                                                                             r_clad,
                                                                             l_mod,
                                                                             num_points_fuel,
                                                                             num_points_ifba,
                                                                             num_points_gap,
                                                                             num_points_clad,
                                                                             num_points_mod)

    num_points, points = get_points_full_pincell(r_fuel,
                                                 r_ifba,
                                                 r_gap,
                                                 r_clad,
                                                 l_mod,
                                                 num_points_fuel,
                                                 num_points_ifba,
                                                 num_points_gap,
                                                 num_points_clad,
                                                 num_points_mod)
    
    node = xml_points(2, # dimension
                      num_points,
                      points)
    spatial_node = node.find("spatial_discretization")
    
    et.ElementTree(node).write(output_path,
                               pretty_print=True,
                               xml_declaration=True)

    if False:
        plt.figure()
        plt.scatter(points[:,0], points[:,1], s=2)
        plt.axes().set_aspect('equal')
        plt.show()
    return

if __name__ == '__main__':
    if (len(sys.argv) != 11):
        print("full_pincell_functions.py [r_fuel r_ifba r_gap r_clad l_mod num_points_fuel num_points_ifba num_points_gap num_points_clad num_points_mod]")
        sys.exit()
    r_fuel = float(sys.argv[1])
    r_ifba = float(sys.argv[2])
    r_gap = float(sys.argv[3])
    r_clad = float(sys.argv[4])
    l_mod = float(sys.argv[5])
    num_points_fuel = int(sys.argv[6])
    num_points_ifba = int(sys.argv[7])
    num_points_gap = int(sys.argv[8])
    num_points_clad = int(sys.argv[9])
    num_points_mod = int(sys.argv[10])
    
    output_full_pincell(r_fuel,
                        r_ifba,
                        r_gap,
                        r_clad,
                        l_mod,
                        num_points_fuel,
                        num_points_ifba,
                        num_points_gap,
                        num_points_clad,
                        num_points_mod)
    
    
                               
    

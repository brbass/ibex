from mesh_functions import *

def get_points_pincell_clad(r_fuel,
                            r_clad,
                            l_mod,
                            num_points_fuel,
                            num_points_clad,
                            num_points_mod):
    mod_points = get_cartesian_mesh_cutout(r_clad,
                                           l_mod,
                                           num_points_mod,
                                           0.7*(r_clad - r_fuel) / num_points_clad)
    fuel_points = get_cylindrical_points(r_fuel,
                                         num_points_fuel,
                                         False)
    clad_points = get_annular_points(r_fuel,
                                     r_clad,
                                     num_points_clad,
                                     True)

    points = np.concatenate((mod_points, fuel_points, clad_points), axis=0)
    num_points = len(points)

    return len(fuel_points), len(clad_points), len(mod_points), num_points, points

def output_pincell_clad(r_fuel,
                        r_clad,
                        l_mod,
                        num_points_fuel,
                        num_points_clad,
                        num_points_mod):
    output_path = "in/pincell_clad_{}_{}_{}_{}_{}_{}.xml".format(r_fuel,
                                                                 r_clad,
                                                                 l_mod,
                                                                 num_points_fuel,
                                                                 num_points_clad,
                                                                 num_points_mod)

    num_total_points_fuel, num_total_points_clad, num_total_points_mod,\
        num_points, points = get_points_pincell_clad(r_fuel,
                                                     r_clad,
                                                     l_mod,
                                                     num_points_fuel,
                                                     num_points_clad,
                                                     num_points_mod)
    
    node = xml_points(2, # dimension
                      num_points,
                      points)
    spatial_node = node.find("spatial_discretization")
    et.SubElement(spatial_node, "number_of_points_fuel").text = str(num_total_points_fuel)
    et.SubElement(spatial_node, "number_of_points_clad").text = str(num_total_points_clad)
    et.SubElement(spatial_node, "number_of_points_moderator").text = str(num_total_points_mod)

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
    if (len(sys.argv) != 7):
        print("pincell_clad_functions.py [r_fuel r_clad l_mod num_points_fuel num_points_clad num_points_mod]")
        sys.exit()
    r_fuel = float(sys.argv[1])
    r_clad = float(sys.argv[2])
    l_mod = float(sys.argv[3])
    num_points_fuel = int(sys.argv[4])
    num_points_clad = int(sys.argv[5])
    num_points_mod = int(sys.argv[6])
    
    output_pincell_clad(r_fuel,
                        r_clad,
                        l_mod,
                        num_points_fuel,
                        num_points_clad,
                        num_points_mod)
    # for i in range(2, 8):
    #     output_pincell_clad(0.4095,
    #                         0.475,
    #                         1.254,
    #                         3 * i,
    #                         i,
    #                         8 * i)
    
    
                               
    

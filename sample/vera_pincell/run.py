from ibex_io import full_run_from_template
from ibex_io import full_save_from_template
from ibex_io import run_multiprocessing_command_no_output
import copy
import sys
import numpy as np
import xml.etree.ElementTree as et
import multiprocessing

def save_point_cases(point_cases, num_procs):
    executable = "python"
    commands = []
    for case in point_cases:
        commands.append([executable, "make_mesh.py --num_radii {} --max_delta {} --initial_delta {} --randomize_start".format(*case)])

    with multiprocessing.Pool(processes=num_procs) as pool:
        pool.starmap(run_multiprocessing_command_no_output, commands)
        pool.close()
        pool.join()

def get_data(case,
             num_procs):
    # Set up data list
    data = {}
    data["executable"] = "ibex"
    data["num_procs"] = num_procs # Number of cases to run simultaneously
    data["parameters"] = ["(POINTS_FILE)",
                          "(TAU)",
                          "(NUM_NEIGHBORS)",
                          "(WEIGHTING)",
                          "(INT_CELLS)",
                          "(RULE)"]
    data["values"] = []
    data["descriptions"] = []
    point_files = []
    descriptions = []
    
    # Set these parameters
    tau_cases = [1.5]
    neighbor_cases = [8]
    weighting_cases = ["full"]
    rules = [5, 4, 3, 2, 1]
    if case == "1b":
        point_cases = [[8, 0.4, 0.2],
                       [8, 0.3, 0.15],
                       [8, 0.2, 0.1],
                       [13, 0.1, 0.05],
                       [13, 0.15, 0.075],
                       [15, 0.05, 0.025],
                       [18, 0.0375, 0.01875],
                       [21, 0.025, 0.0125],
                       [24, 0.01875, 0.009375],
                       [27, 0.0125, 0.00625],
                       [25, 0.01, 0.001]]
    else: # case == "1e"
        point_cases = [[18, 0.06, 0.006],
                       [19, 0.04, 0.004],
                       [20, 0.03, 0.003],
                       [22, 0.02, 0.002],
                       [24, 0.015, 0.0015],
                       [25, 0.01, 0.001]]
    save_point_cases(point_cases, num_procs)
    for dat in point_cases:
        point_files.append("vera_mesh_{}_{}_{}_0.4101_1.26_rand.xml".format(*dat))
        descriptions.append("mult{}-{}-{}-r".format(*dat))
        
    # Get point files and data
    points = []
    int_cells = []
    for point_file in point_files:
        num_points = int(et.parse(point_file).getroot().find("spatial_discretization").findtext("number_of_points"))
        points.append(num_points)
        num_cells = int(4 * np.sqrt(num_points))
        int_cells.append(num_cells)
        
    # Get permutations of parameters
    for rule in rules:
        for point_file, description, num_points, num_cells in zip(point_files, descriptions, points, int_cells):
            for tau in tau_cases:
                for neighbors in neighbor_cases:
                    for weighting in weighting_cases:
                        mem = 0.0006 * num_points * np.power(4., rule - 3)
                        if mem < 100: # Set this to the maximum memory in GB
                            data["values"].append([point_file,
                                                   tau,
                                                   neighbors,
                                                   weighting,
                                                   num_cells,
                                                   rule])
                            data["descriptions"].append([description,
                                                         tau,
                                                         neighbors,
                                                         weighting,
                                                         num_cells,
                                                         rule])
    
    data["prefix"] = "test"
    data["postfix"] = ".xml"
    data["template_filename"] = "template_{}.xml".format(case)
    
    return data

def run(case,
        num_procs,
        run = False):
    # Get data
    data = get_data(case, num_procs)
    
    # Run cases
    if run:
        input_filenames = full_run_from_template(data,
                                                 True) # Save input files
    else:
        input_filenames = full_save_from_template(data,
                                                  True) # Save input files
if __name__ == '__main__':
    if sys.version_info[0] < 3:
        print("This script requires Python 3")
        sys.exit(1)

    run_problems = True
    num_cases = 4
    for case in ["1b", "1e"]:
        run(case, num_cases, run_problems)

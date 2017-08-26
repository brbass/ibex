import sys
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import tri
import xml.etree.ElementTree as et
import itertools
import multiprocessing
import subprocess
import os

# Unpack an indexed array
def unpack2(data, num1, num2):
    if (len(data) != num1 * num2):
        print("error: data wrong size")
    newdata = np.zeros((num1, num2))
    for i in range(num1):
        for j in range(num2):
            newdata[i, j] = data[j + num2 * i]
    return newdata
def unpack3(data, num1, num2, num3):
    if (len(data) != num1 * num2 * num3):
        print("error: data wrong size")
    newdata = np.zeros((num1, num2, num3))
    for i in range(num1):
        for j in range(num2):
            for k in range(num3):
                newdata[i, j, k] = data[k + num3 * (j + num2 * i)]
    return newdata

# Get scalar data from XML
def get_scalar(data_name, dtype, node):
    return dtype(node.findtext(data_name))

# Get a 1D vector from XML
def get_vector(data_name, dtype, node):
    return np.fromstring(node.findtext(data_name), dtype=dtype, sep="\t")

# Put pertinant output data into a struct
def get_data(file_path):
    # Read in output file
    output_node = et.parse(file_path).getroot()
    
    # Initialize data struct
    data = {}

    # Energy discretization
    child_node = output_node.find("energy_discretization")
    data["number_of_groups"] = get_scalar("number_of_groups", int, child_node)
    num_groups = data["number_of_groups"]
    
    # Angular discretization
    child_node = output_node.find("angular_discretization")
    data["number_of_moments"] = get_scalar("number_of_moments", int, child_node)
    num_moments = data["number_of_moments"]
    data["number_of_ordinates"] = get_scalar("number_of_ordinates", int, child_node)
    num_ordinates = data["number_of_ordinates"]
    
    # Spatial discretization
    child_node = output_node.find("spatial_discretization")
    data["dimension"] = get_scalar("dimension", int, child_node)
    dimension = data["dimension"]
    data["number_of_points"] = get_scalar("number_of_points", int, child_node)
    num_points = data["number_of_points"]
    data["points"] = unpack2(get_vector("points", float, child_node),
                             num_points,
                             dimension)
    
    # Spatial discretization options
    child_node = child_node.find("options")
    data["integration_ordinates"] = get_scalar("integration_ordinates", int, child_node)
    data["supg"] = get_scalar("include_supg", bool, child_node)
    data["weighting"] = get_scalar("weighting", str, child_node)
    data["tau_scaling"] = get_scalar("tau_scaling", str, child_node)

    # Solver
    child_node = output_node.find("solver")
    data["total_iterations"] = get_scalar("total_iterations", int, child_node)
    data["coefficients"] = unpack3(get_vector("coefficients", float, child_node),
                                   num_points,
                                   num_moments,
                                   num_groups)
    data["k_eigenvalue"] = get_scalar("k_eigenvalue", float, child_node)

    # Solver values
    child_node = child_node.find("values")
    data["phi"] = unpack3(get_vector("phi", float, child_node),
                          num_points,
                          num_moments,
                          num_groups)

    # Timing
    timing = {}
    child_node = output_node.find("timing")
    timing["spatial_initialization"] = get_scalar("spatial_initialization", float, child_node)
    timing["sweep_initialization"] = get_scalar("sweep_initialization", float, child_node)
    timing["solve"] = get_scalar("solve", float, child_node)
    data["timing"] = timing
    
    return data

# Data should contain parameters, a list of strings to be replaced,
# and values, a list for each of these parameters.
# Ex: parameters = ["num_points", "num_groups"], values = [[120, 240], [2, 3]]
# Assumes values are scalars or strings
def save_from_template(data):
    # Get template string
    template_filename = data["template_filename"]
    template_file = open(template_filename, "r")
    template_string = template_file.read()
    template_file.close()
    
    # Get data
    parameters = data["parameters"]
    num_parameters = len(parameters)
    values = data["values"]
    prefix = data["prefix"]
    postfix = data["postfix"]
    
    # Get all permutations of data
    perms = list(itertools.product(*values))

    # Get input files
    input_filenames = []
    for perm in perms:
        # Get input filename and string
        input_string = template_string
        input_filename = prefix
        for parameter, value in zip(parameters, perm):
            input_string = input_string.replace(parameter, str(value))
            input_filename += "_{}".format(value)
        input_filename += postfix
        input_filenames.append(input_filename)

        # Save input file
        input_file = open(input_filename, "w")
        input_file.write(input_string)
        input_file.close()
    
    return input_filenames

def run_multiprocessing_command(executable,
                                arguments):
    # Change to working directory
    # starting_directory = os.getcwd()
    # os.chdir(directory)
    
    # Get command
    pid = multiprocessing.current_process().name
    command = "{} {}".format(executable, arguments)

    # Run command
    print("start \"{}\" on process {}".format(command, pid))
    subprocess.call([command], shell=True)
    print("end \"{}\" on process {}".format(command, pid))

    # Switch back to starting directory
    # os.chdir(starting_directory)

def run_from_template(data):
    # Get commands
    executable = data["executable"]
    # directory = data["directory"]
    input_filenames = save_from_template(data)
    input_commands = [[executable, input_filename] for input_filename in input_filenames]
    
    # Create multiprocessing pool
    num_procs = data["num_procs"]
    pool = multiprocessing.Pool(processes=num_procs)
    
    # Run problems
    pool.starmap(run_multiprocessing_command, input_commands)
    

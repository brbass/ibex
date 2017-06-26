import sys
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import tri
import xml.etree.ElementTree as et

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
    
    

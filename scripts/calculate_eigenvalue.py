import sys
import numpy as np
import numpy.linalg as npl
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

# Get scalar data from XML
def get_scalar(data_name, dtype, node):
    return dtype(node.findtext(data_name))

# Get a 1D vector from parent XML node
def get_vector(data_name, dtype, node):
    return np.fromstring(node.findtext(data_name), dtype=dtype, sep="\t")

def output_eigenvalues(input_file):
    # Get file as XML
    node = et.parse(input_file).getroot()
    material_nodes = node.find("materials")

    # Get number of energy groups
    num_groups = get_scalar("number_of_groups", int, node.find("energy_discretization"))
    
    # Load each material in turn
    for material_node in material_nodes.iter("material"):
        sigma_t = get_vector("sigma_t", float, material_node)
        sigma_s = unpack2(get_vector("sigma_s", float, material_node),
                          num_groups,
                          num_groups)
        sigma_f = get_vector("sigma_f", float, material_node)
        nu = get_vector("nu", float, material_node)
        chi = get_vector("chi", float, material_node)
        
        mat1 = np.zeros((num_groups, num_groups))
        mat2 = np.zeros((num_groups, num_groups))
        for g1 in range(num_groups):
            for g2 in range(num_groups):
                if g1 == g2:
                    mat1[g1, g2] += sigma_t[g1]
                mat1[g1, g2] -= sigma_s[g1, g2]
                mat2[g1, g2] = chi[g1] * nu[g2] * sigma_f[g2]
        eig_mat = np.matmul(npl.inv(mat1), mat2)
        eig = npl.eig(eig_mat)
        print(eig[0])

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("calculate_eigenvalue.py [ibex_input_file.xml]")
        exit()
    output_eigenvalues(sys.argv[1])

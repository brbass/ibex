import os, sys, subprocess, itertools, multiprocessing, functools
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm as cm
from ibex_io import get_data

def plot_flux(file_path,
              file_out = ""):
    data = get_data(file_path)
    dimension = data["dimension"]
    if dimension == 1:
        plot_flux_1d(data,
                     file_out)
    elif dimension == 2:
        plot_flux_2d(data,
                     file_out)

def plot_flux_1d(data,
                 file_out):
    num_groups = data["number_of_groups"]
    num_points = data["number_of_points"]
    num_moments = data["number_of_moments"]
    phi = data["phi"]
    points = data["points"]

    plt.figure()
    for g in range(num_groups):
        plt.plot(points[:, 0], phi[:, 0, g], label="g{}".format(g))
    plt.xlabel("x")
    plt.ylabel("$\phi(x)$")
    plt.legend()
    plt.tight_layout()
    
    if (file_out == ""):
        plt.show()
    else:
        plt.savefig(file_out + "_flux_1d.pdf", bbox_inches='tight')
        plt.close()

def plot_flux_2d(data,
                 file_out):
    num_groups = data["number_of_groups"]
    num_points = data["number_of_points"]
    num_moments = data["number_of_moments"]
    phi = data["phi"]
    points = data["points"]
    
    for g in range(num_groups):
        plt.figure(g)
        plt.tripcolor(points[:, 0], points[:, 1], phi[:, 0, g], shading='gouraud')
        plt.xlabel("x")
        plt.ylabel("y")
        plt.title("group " + str(g))
        plt.colorbar(label=r"$\phi(x, y)$")
        plt.axis('equal')
        plt.tight_layout()
        
    if (file_out == ""):
        plt.show()
    else:
        for g in range(num_groups):
            plt.figure(g)
            plt.savefig(file_out + "_flux_2d_" + str(g) + ".pdf", bbox_inches='tight')
            plt.close()

if __name__ == '__main__':
    if (len(sys.argv) != 2 and len(sys.argv) != 3):
        print("ibex_plot.py [file_in (file_out)]")
        sys.exit()
    file_path = sys.argv[1]
    if (len(sys.argv) == 3):
        file_out = sys.argv[2]
        plot_flux(file_path,
                  file_out)
    else:
        plot_flux(file_path)
    
            

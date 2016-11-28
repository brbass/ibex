import numpy as np
import sys
from matplotlib import pyplot as plt
import scipy.linalg as spl
import scipy.integrate as spi
import scipy.signal as sps
from basis import *

def compare_strong_weak():
    def func(x):
        return sps.sawtooth(2 * 2 * np.pi * x, width=0.5) + 1 # sawtooth
        # return sps.square(2 * np.pi * ((2 * x) - 0.25)) # square
        # return np.sin(2 * np.pi * (3.5 * x))
    
    num_points = 17
    points = np.linspace(0, 1, num=num_points)
    # basis = Compact_Gaussian(1.0,
    #                          points)
    # weight = Compact_Gaussian(1.0,
    #                           points)
    basis = Linear_MLS(2,
                       points)
    weight = Linear_MLS(2,
                        points)
    num_plot_points = (num_points-1)*20+1
    plot_points = np.linspace(0, 1, num=num_plot_points)
    ana_plot = func(plot_points)

    # Strong interpolation
    
    mat = np.zeros((num_points, num_points))
    for i in range(num_points):
        for j, point in enumerate(points):
            mat[i, j] = basis.val(i, point)
    rhs = func(points)
    strong_sol = spl.solve(mat, rhs)

    strong_plot = np.zeros(plot_points.size)
    for i, point in enumerate(plot_points):
        for j in range(num_points):
            strong_plot[i] += strong_sol[j] * basis.val(j, point)

    # Weak interpolation
    
    mat = np.zeros((num_points, num_points))
    for i in range(num_points):
        for j in range(num_points):
            nonzero, limits = get_limits(i,
                                         j,
                                         basis,
                                         weight)
            if nonzero:
                def integrand(x):
                    bas = basis.val(i, x)
                    wei = weight.val(j, x)
                    return bas * wei
                mat[j, i], err = spi.quad(integrand, limits[0], limits[1])
        limits = weight.limits(i)
        def integrand(x):
            wei = weight.val(i, x)
            f = func(x)
            return wei * f
        rhs[i], err = spi.quad(integrand, limits[0], limits[1])
    weak_sol = spl.solve(mat, rhs)
    weak_plot = np.zeros(plot_points.size)
    for i, point in enumerate(plot_points):
        for j in range(num_points):
            weak_plot[i] += weak_sol[j] * basis.val(j, point)

    # Plot comparison
    
    plt.figure()
    plt.plot(plot_points, strong_plot, label="strong")
    plt.plot(plot_points, weak_plot, label="weak")
    plt.plot(plot_points, ana_plot, label="analytic")
    plt.legend()

    # Plot constituent parts
    
    plt.figure()
    for i in range(num_points):
        limits = basis.limits(i)
        points_temp = np.linspace(0, 1, num=num_plot_points)
        part_sol = np.zeros(num_plot_points)
        for j, point in enumerate(points_temp):
            part_sol[j] = strong_sol[i] * basis.val(i, point)
        plt.plot(plot_points, part_sol)

    # Plot plain basis functions
    
    plt.figure()
    for i in range(num_points):
        limits = basis.limits(i)
        points_temp = np.linspace(0, 1, num=num_plot_points)
        part_sol = np.zeros(num_plot_points)
        for j, point in enumerate(points_temp):
            part_sol[j] = basis.val(i, point)
        plt.plot(plot_points, part_sol)

    plt.show()
    
if __name__ == '__main__':
    compare_strong_weak()



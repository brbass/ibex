import numpy as np
import sys
from matplotlib import pyplot as plt
import scipy.linalg as spl
import scipy.integrate as spi
import scipy.signal as sps
from basis import *

colors = ['#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854']

def compare_strong_weak(num_points, func_type, basis_type, shape):
    if func_type == "saw":
        def func(x):
            return sps.sawtooth(2 * 2 * np.pi * x, width=0.5) + 1 # sawtooth
    elif func_type == "square":
        def func(x):
            return sps.square(2 * np.pi * ((2 * x) - 0.25)) + 1 # square
    else:
        def func(x):
            return np.sin(2 * np.pi * (3.5 * x))
        
    description = "../figs/{}_{}_{}_{}".format(func_type,
                                               basis_type,
                                               num_points,
                                               shape)
    points = np.linspace(0, 1, num=num_points)

    if basis_type == "gauss":
        basis = Compact_Gaussian(shape,
                                 points)
        weight = Compact_Gaussian(shape,
                                  points)
    elif basis_type == "mq":
        basis = Multiquadric(shape,
                             points)
        weight = Multiquadric(shape,
                              points)
    elif basis_type == "mls":
        basis = Linear_MLS(shape,
                           points)
        weight = Linear_MLS(shape,
                            points)
    num_plot_points = (num_points-1)*20+1
    plot_points = np.linspace(0, 1, num=num_plot_points)
    ana_plot = func(plot_points)

    # Strong interpolation
    
    mat = np.zeros((num_points, num_points))
    for i in range(num_points):
        for j, point in enumerate(points):
            mat[j, i] = basis.val(i, point)
    rhs = func(points)
    strong_sol = spl.solve(mat, rhs)
    const_sol = spl.solve(mat, np.ones(num_points))
    lin_sol = spl.solve(mat, np.linspace(0, 1, num=num_points))
    
    strong_plot = np.zeros(plot_points.size)
    const_plot = np.zeros(plot_points.size)
    lin_plot = np.zeros(plot_points.size)
    for i, point in enumerate(plot_points):
        for j in range(num_points):
            strong_plot[i] += strong_sol[j] * basis.val(j, point)
            const_plot[i] += const_sol[j] * basis.val(j, point)
            lin_plot[i] += lin_sol[j] * basis.val(j, point)
    
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
    
    fig, ax1 = plt.subplots()
    # ax2 = ax1.twinx()
    ax1.set_xlabel(r"$x$")
    ax1.set_ylabel(r"$f(x)$")
    # ax2.set_ylabel(r"$f(x) - \tilde{f}(x)$")
    ln1 = ax1.plot(plot_points, strong_plot, label="strong", color=colors[0])
    ln2 = ax1.plot(plot_points, weak_plot, label="weak", color=colors[1])
    ln3 = ax1.plot(plot_points, ana_plot, label="analytic", color=colors[2])
    # ln4 = ax2.plot(plot_points, strong_plot - ana_plot, label="err_strong", color=colors[3])
    # ln5 = ax2.plot(plot_points, weak_plot - ana_plot, label="err_weak", color=colors[4])
    ax1.legend(loc='lower left')
    # ax2.legend(loc='lower right')
    ax1.grid()
    plt.savefig(description + "_interp.pdf", bbox_inches='tight')
    plt.close()
    
    # Plot constituent parts
    
    plt.figure()
    for i in range(num_points):
        limits = basis.limits(i)
        points_temp = np.linspace(0, 1, num=num_plot_points)
        part_sol = np.zeros(num_plot_points)
        for j, point in enumerate(points_temp):
            part_sol[j] = strong_sol[i] * basis.val(i, point)
        plt.plot(plot_points, part_sol, color='k')
    plt.xlabel(r"$x$")
    plt.ylabel(r"$\Gamma_i(x)a_i$")
    plt.axes().grid()
    plt.savefig(description + "_partial.pdf", bbox_inches='tight')
    plt.close()

    # Plot plain basis functions
    
    plt.figure()
    for i in range(num_points):
        limits = basis.limits(i)
        points_temp = np.linspace(0, 1, num=num_plot_points)
        part_sol = np.zeros(num_plot_points)
        for j, point in enumerate(points_temp):
            part_sol[j] = basis.val(i, point)
        plt.plot(plot_points, part_sol, color='k')
    plt.xlabel(r"$x$")
    plt.ylabel(r"$\Gamma_i(x)$")
    plt.axes().grid()
    plt.savefig(description + "_basis.pdf", bbox_inches='tight')
    plt.close()

    # Plot constant and linear solution
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.set_xlabel(r"$x$")
    ax1.set_ylabel(r"$f(x)$")
    ax2.set_ylabel(r"$err(f(x))$")
    ln1 = ax1.plot(plot_points, const_plot, label="const", color=colors[0])
    ln2 = ax1.plot(plot_points, lin_plot, label="linear", color=colors[1])
    ln4 = ax2.semilogy(plot_points, np.abs(const_plot - 1), label="err_const", color=colors[2])
    ln5 = ax2.semilogy(plot_points, np.abs(lin_plot - np.linspace(0, 1, num=num_plot_points)), label="err_linear", color=colors[3])
    ax1.legend(loc='lower center')
    ax2.legend(loc='lower right')
    ax1.grid()
    plt.savefig(description + "_const.pdf", bbox_inches='tight')
    plt.close()
    
if __name__ == '__main__':
    compare_strong_weak(17, "square", "gauss", 1.0)
    compare_strong_weak(17, "square", "mls", 2)



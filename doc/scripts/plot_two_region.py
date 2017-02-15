import numpy as np
from matplotlib import pyplot as plt

def plot_two_region(number_of_points,
                    positions,
                    radii,
                    point_types,
                    filename_out):
    include_radii = True
    include_centers = False
    bound_color='k'
    basis_color='b'
    weight_color='r'
    point_color='g'
    bound_width = 2.
    num_plot = 50
    theta = np.linspace(0, 2*np.pi, num_plot)
    bound_min = np.zeros(num_plot)
    bound_max = np.ones(num_plot)
    bound_mid = 0.5 * np.ones(num_plot)
    bound_range = np.linspace(0, 1, num_plot)
    
    # Create figure
    plt.figure()
    plt.axis('equal')
    plt.xlim(-1, 2)
    plt.ylim(-1, 2)
    plt.xlabel("x")
    plt.ylabel("y")
    
    # Plot boundaries
    plt.plot(bound_min, bound_range, color=bound_color, linewidth=bound_width, label="boundaries")
    plt.plot(bound_mid, bound_range, color=bound_color, linewidth=bound_width)
    plt.plot(bound_max, bound_range, color=bound_color, linewidth=bound_width)
    plt.plot(bound_range, bound_min, color=bound_color, linewidth=bound_width)
    plt.plot(bound_range, bound_max, color=bound_color, linewidth=bound_width)

    # Plot points
    if include_centers:
        plt.scatter(positions[:, 0], positions[:, 1], color=point_color, marker='o', label="centers")
    
    # Plot radii
    if include_radii:
        have_included_basis_label = False
        have_included_weight_label = False
        for i in range(number_of_points):
            position = positions[i]
            radius = radii[i]
            
            x_vals = position[0] + radius * np.cos(theta)
            y_vals = position[1] + radius * np.sin(theta)
            if point_types[i] == 0:
                if have_included_basis_label:
                    plt.plot(x_vals, y_vals, color=basis_color)
                else:
                    plt.plot(x_vals, y_vals, color=basis_color, label="basis")
                    have_included_basis_label = True
            else:
                if have_included_weight_label:
                    plt.plot(x_vals, y_vals, color=weight_color)
                else:
                    plt.plot(x_vals, y_vals, color=weight_color, label="weight")
                    have_included_weight_label = True

    # Save figure
    plt.legend()
    plt.savefig(filename_out, bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    if False:
        number_of_points = 6
        positions = np.array([[0.5, 0.5], # weight function
                              [0.6, 0.6], # basis contained in weight
                              [0.45, 0.55], # weight contained in basis
                              [0.35, 0.3], # intersection after centerline of basis
                              [0.2, 0.7], # intersection before centerline of basis
                              [0.7, 0.2]]) # basis includes boundary
        radii = [0.3,
                 0.1,
                 0.4,
                 0.1,
                 0.1,
                 0.55]
        
        point_types = [1, 0, 0, 0, 0, 0]
        plot_two_region(number_of_points,
                        positions,
                        radii,
                        point_types,
                        "../figs/weight_test_1.pdf")
    if False:
        number_of_points = 6
        positions = np.array([[0., 0.], 
                              [0.25, 0.3], 
                              [0.2, 0.5], 
                              [0.3, 0.2], 
                              [0.4, 0.], 
                              [0.1, 0.15]]) 
        radii = [0.4,
                 0.5,
                 0.3,
                 0.1,
                 0.15,
                 0.05]
        
        point_types = [1, 0, 0, 0, 0, 0]
        plot_two_region(number_of_points,
                        positions,
                        radii,
                        point_types,
                        "../figs/weight_test_2.pdf")

    if False:
        number_of_points = 6
        positions = np.array([[0.9, 0.5], 
                              [0.9, 0.9], 
                              [1., 0.4], 
                              [0.7, 0.2], 
                              [0.8, 0.0], 
                              [0.2, 0.6]]) 
        radii = [0.4,
                 0.2,
                 0.25,
                 0.15,
                 0.2,
                 0.4]
        
        point_types = [1, 0, 0, 0, 0, 0]
        plot_two_region(number_of_points,
                        positions,
                        radii,
                        point_types,
                        "../figs/weight_test_3.pdf")
    if True:
        number_of_points = 2
        positions = np.array([[-1./3., -1./3.],
                              [1./9., -1]])
        radii = np.array([4./9.,
                          0.4969039949999532])
        point_types = [0, 1]
        plot_two_region(number_of_points,
                        positions,
                        radii,
                        point_types,
                        "./test.pdf")

import numpy as np
from matplotlib import pyplot as plt

colors = ['#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854']

def plot_lens(d, r1, r2):
    num_vals = 50
    phivals = np.linspace(0, 2*np.pi, num=num_vals)
    xinter = (r1*r1 - r2*r2 + d*d) / (2*d)
    fig, ax = plt.subplots()
    for i in range(2):
        if i == 0:
            xvals = r1 * np.cos(phivals)
            yvals = r1 * np.sin(phivals)
            xfvals = np.linspace(d-r2, xinter, num=num_vals)
            yfvals1 = np.sqrt(r2 * r2 - np.power(xfvals - d, 2))
        else:
            xvals = d + r2 * np.cos(phivals)
            yvals = r2 * np.sin(phivals)
            xfvals = np.linspace(xinter, r1, num=num_vals)
            yfvals1 = np.sqrt(r1 * r1 - np.power(xfvals, 2))
        ax.plot(xvals, yvals, color=colors[i])
        yfvals2 = -yfvals1
        ax.fill_between(xfvals, yfvals2, yfvals1, color=colors[3])
    ax.set_aspect(1)
    ax.grid()
    plt.savefig("../figs/lens_{}_{}_{}.pdf".format(d, r1, r2),
                bbox_inches='tight')
    
if __name__ == '__main__':
    plot_lens(2, 1, 2.5)
    plot_lens(2.5, 1, 2.5)
    plot_lens(2, 1, 3)
        

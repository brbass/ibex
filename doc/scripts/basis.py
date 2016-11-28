import numpy as np
import sys
from matplotlib import pyplot as plt

# Basis and weight functions

def get_limits(i,
               j,
               basis,
               weight):
    baslim = basis.limits(i)
    weilim = weight.limits(j)
    newlim = np.zeros(2)
    
    if baslim[1] < weilim[0] or baslim[0] > weilim[1]:
        return False, newlim
    newlim[0] = np.maximum(baslim[0], weilim[0])
    newlim[1] = np.minimum(baslim[1], weilim[1])

    return True, newlim

class RBF:
    def __init__(self,
                 shape,
                 points):
        self.shape = shape / (points[1] - points[0])
        self.points = points
        self.compact = False
    def limits(self,
               i):
        return [self.points[0], self.points[-1]]

class Compact_RBF(RBF):
    def __init__(self,
                 shape,
                 points,
                 max_distance):
        RBF.__init__(self,
                     shape,
                     points)
        self.max_distance = max_distance
        self.compact = True
        self.limit = np.zeros((len(self.points), 2), dtype=float)
        for i, point in enumerate(self.points):
            self.limit[i, 0] = np.amax([0., self.points[i] - self.max_distance / self.shape])
            self.limit[i, 1] = np.amin([self.points[-1], self.points[i] + self.max_distance / self.shape])
    def limits(self,
               i):
        return self.limit[i, :]
        
class Multiquadric(RBF):
    def __init__(self,
                 shape,
                 points):
        RBF.__init__(self,
                     shape,
                     points)

    def val(self,
            i,
            x):
        d = x - self.points[i]
        return np.sqrt(1 + np.power(d * self.shape, 2))
    def dval(self,
             i,
             x):
        d = x - self.points[i]
        return np.power(self.shape, 2) * d / np.sqrt(1 + np.power(d * self.shape, 2))

class Compact_Gaussian(Compact_RBF):
    def __init__(self,
                 shape,
                 points):
        Compact_RBF.__init__(self,
                             shape,
                             points,
                             3.)
    
    def val(self,
            i,
            x):
        d = x - self.points[i]
        if np.abs(d * self.shape) < self.max_distance:
            return np.exp(-np.power(d * self.shape, 2))
        else:
            return 0.
        
    def dval(self,
             i,
             x):
        d = x - self.points[i]
        if np.abs(d * self.shape) < self.max_distance:
            return -2 * np.power(self.shape, 2) * d * np.exp(-np.power(d * self.shape, 2))
        else:
            return 0.
    def ddval(self,
              i,
              x):
        d = x - self.points[i]
        
        if np.abs(d * self.shape) < self.max_distance:
            return 2. * np.exp(2) * np.power(self.shape, 2) * (2. * np.exp(2) * np.power(d, 2) - 1) * np.exp(-np.power(d * self.shape, 2))
        else:
            return 0.
        
class Gaussian(RBF):
    def __init__(self,
                 shape,
                 points):
        RBF.__init__(self,
                     shape,
                     points)
    def val(self,
            i,
            x):
        d = x - self.points[i]
        return np.exp(-np.power(d * self.shape, 2))
    def dval(self,
             i,
             x):
        d = x - self.points[i]
        return -2 * np.power(self.shape, 2) * d * np.exp(-np.power(d * self.shape, 2))

class Wendland(Compact_RBF):
    def __init__(self,
                 shape,
                 points):
        Compact_RBF.__init__(self,
                             shape,
                             points,
                             1.)
        
    def val(self,
            i,
            x):
        d = np.abs(x - self.points[i]) * self.shape
        if d < 1:
            return np.power(1 - d, 4) * (1 + 4 * d)
        else:
            return 0.
        
    def dval(self,
             i,
             x):
        d = np.abs(x - self.points[i]) * self.shape
        if d < 1:
            return 20 * self.shape * d * np.power(d - 1, 3)
        else:
            return 0.

class Linear_MLS(Compact_RBF):
    def __init__(self,
                 num_other_points,
                 points):
        self.num_polynomials = 2
        self.num_other_points = num_other_points
        self.points = points
        self.num_points = len(points)
        self.dx = points[1] - points[0]
        self.bandwidth = 1. / (self.dx * (self.num_other_points + 0.1))
        self.shape = self.bandwidth
        Compact_RBF.__init__(self,
                             self.dx,
                             self.points,
                             1. / self.bandwidth)
    def weight(self,
               i,
               x):
        d = self.bandwidth * (x - self.points[i])
        d1 = np.abs(d)
        if d1 <= 1:
            return np.power(1. - d1, 3) * (1. + 3. * d1)
        else:
            return 0.
        
    def dweight(self,
                i,
                x):
        d = self.bandwidth * (x - self.points[i])
        d1 = np.abs(d)
        if d1 <= 1:
            return -12 * np.sign(d) * self.bandwidth * np.power(d1 - 1, 2) * d1
        else:
            return 0.
        
    def val(self,
            i,
            x):
        xi = self.points[i]
        if np.abs(self.bandwidth * (x - xi)) > 1:
            return 0.
        polyval = np.array([1, x])
        a = np.zeros(3)
        nearest_point = int(np.round(x / self.dx))
        first_point = nearest_point - self.num_other_points - 1
        last_point = nearest_point + self.num_other_points + 1
        if first_point < 0:
            first_point = 0
        if last_point >= self.num_points:
            last_point = self.num_points - 1
        for j in range(first_point, last_point + 1):
            xj = self.points[j]
            w = self.weight(j, x)
            a[0] += w
            a[1] += w * xj
            a[2] += w * xj * xj
        denom = a[0]*a[2] - a[1]*a[1]
        ainv = np.array([[a[2]/denom, -a[1]/denom],
                         [-a[1]/denom, a[0]/denom]])
        b = self.weight(i, x) * np.array([1, xi])
        return np.dot(polyval, np.dot(ainv, b))
    
    def dval(self,
             i,
             x):
        xi = self.points[i]
        if np.abs(self.bandwidth * (x - xi)) > 1:
            return 0.
        polyval = np.array([1, x])
        dpolyval = np.array([0, 1])
        a = np.zeros(3)
        da = np.zeros((2,2))
        nearest_point = int(np.round(x / self.dx))
        first_point = nearest_point - self.num_other_points - 1
        last_point = nearest_point + self.num_other_points + 1
        if first_point < 0:
            first_point = 0
        if last_point >= self.num_points:
            last_point = self.num_points - 1
        for j in range(first_point, last_point + 1):
            xj = self.points[j]
            w = self.weight(j, x)
            a[0] += w
            a[1] += w * xj
            a[2] += w * xj * xj
            dw = self.dweight(j, x)
            da[0, 0] += dw
            da[1, 0] += dw * xj
            da[0, 1] += dw * xj
            da[1, 1] += dw * xj * xj
        denom = a[0]*a[2] - a[1]*a[1]
        ainv = np.array([[a[2]/denom, -a[1]/denom],
                         [-a[1]/denom, a[0]/denom]])
        polypoint = np.array([1, xi])
        b = self.weight(i, x) * polypoint
        db = self.dweight(i, x) * polypoint
        dainv = -np.dot(ainv, np.dot(da, ainv))
        t1 = np.dot(dpolyval, np.dot(ainv, b))
        t2 = np.dot(polyval, np.dot(dainv, b))
        t3 = np.dot(polyval, np.dot(ainv, db))
        return t1 + t2 + t3

class MLS(Compact_RBF):
    def __init__(self,
                 num_polynomials,
                 num_other_points,
                 points):
        self.num_polynomials = num_polynomials
        self.num_other_points = num_other_points
        self.points = points
        self.num_points = len(points)
        self.dx = points[1] - points[0]
        self.bandwidth = 1. / (self.dx * (self.num_other_points + 0.1))
        self.shape = self.bandwidth
        Compact_RBF.__init__(self,
                             self.dx,
                             self.points,
                             1. / self.bandwidth)
    def weightd(self,
                d):
        d1 = np.abs(d)
        if d1 <= 1:
            return np.power(1. - d1, 3) * (1. + 3. * d1)
        else:
            return 0.
    def dweightd(self,
                 d):
        d1 = np.abs(d)
        if d1 <= 1:
            return -12 * np.sign(d) * self.bandwidth * np.power(d1 - 1, 2) * d1
        else:
            return 0.
    def weight(self,
               i,
               x):
        return self.weightd(self.bandwidth * (x - self.points[i]))
    def dweight(self,
                i,
                x):
        return self.dweightd(self.bandwidth * (x - self.points[i]))
    def get_points(self,
                   x):
        nearest_point = int(np.round(x / self.dx))
        local_points = []
        for i in range(nearest_point - self.num_other_points - 1, nearest_point + self.num_other_points + 1):
            if (i >= 0 and i < self.num_points
                and np.abs(x - self.points[i]) < self.bandwidth):
                local_points.append([i])
        return np.array(local_points, dtype=int)
    def amat(self,
             x):
        local_points = self.get_points(x)
        mat = np.zeros((self.num_polynomials, self.num_polynomials))
        
        for i in local_points:
            w = self.weight(i,
                            x)
            for j in range(self.num_polynomials):
                for k in range(self.num_polynomials):
                    mat[j, k] = mat[j, k] + np.power(self.points[i], j + k) * w
        return mat
    def damat(self,
              x):
        local_points = self.get_points(x)
        mat = np.zeros((self.num_polynomials, self.num_polynomials))
        
        for i in local_points:
            dw = self.dweight(i,
                              x)
            for j in range(self.num_polynomials):
                for k in range(self.num_polynomials):
                    mat[j, k] = mat[j, k] + np.power(self.points[i], j + k) * dw
        return mat
    def ainv(self,
             x):
        return np.linalg.inv(self.amat(x))
    def dainv(self,
              x):
        ainvval = self.ainv(x)
        daval = self.damat(x)
        return -1. * np.dot(ainvval, np.dot(daval, ainvval))
    def poly(self,
             x):
        vec = np.zeros(self.num_polynomials)
        for i in range(self.num_polynomials):
            vec[i] = np.power(x, i)
        return vec
    def dpoly(self,
              x):
        vec = np.zeros(self.num_polynomials)
        for i in range(1, self.num_polynomials):
            vec[i] = i * np.power(x, i - 1)
        return vec
    def bvec(self,
             i,
             x):
        return self.poly(self.points[i]) * self.weight(i,
                                                       x)
    def dbvec(self,
              i,
              x):
        return self.poly(self.points[i]) * self.dweight(i,
                                                        x)
        
    def val(self,
            i,
            x):
        polyval = self.poly(x)
        ainvval = self.ainv(x)
        bvecval = self.bvec(i,
                            x)
        return np.dot(polyval, np.dot(ainvval, bvecval))
    def dval(self,
             i,
             x):
        polyval = self.poly(x)
        dpolyval = self.dpoly(x)
        aval = self.amat(x)
        daval = self.damat(x)
        ainvval = np.linalg.inv(aval)
        bval = self.bvec(i,
                         x)
        dbval = self.dbvec(i,
                           x)
        dainvval = -np.dot(ainvval, np.dot(daval, ainvval))
        t1 = np.dot(dpolyval, np.dot(ainvval, bval))
        t2 = np.dot(polyval, np.dot(dainvval, bval))
        t3 = np.dot(polyval, np.dot(ainvval, dbval))
        return t1 + t2 + t3

class Constant(Compact_RBF):
    def __init__(self,
                 shape,
                 points):
        Compact_RBF.__init__(self,
                             shape,
                             points,
                             1.)
        
    def val(self,
            i,
            x):
        d = np.abs(x - self.points[i]) * self.shape
        if d < 1:
            return 1.
        else:
            return 0.
        
    def dval(self,
             i,
             x):
        d = np.abs(x - self.points[i]) * self.shape
        if d < 1:
            return 0.
        else:
            return 0.

if __name__ == '__main__':
    points = np.linspace(0, 1, 10)
    shape = 2.0
    # funcs = [Wendland(shape,
    #                   points),
    #          Multiquadric(shape,
    #                       points),
    #          Gaussian(shape,
    #                   points),
    #          Constant(shape,
    #                   points),
    #          Compact_Gaussian(shape,
    #                           points)]
    # desc = ["wend", "mult", "gauss", "const", "com_gauss", "supg_gauss"]
    funcs = [MLS(2,
                 3,
                 points)]
    desc = ["mls"]
    col = ['#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854', 'k']
    plot_points = np.linspace(0, 1, 200)
    for i, func in enumerate(funcs):
        for j, point in enumerate(points):
            vals = np.array([func.val(j, x) for x in plot_points])
            dvals = np.array([func.dval(j, x) for x in plot_points])
            if j == 0:
                plt.figure(0)
                plt.plot(plot_points, vals, color=col[i], label=desc[i])
                plt.figure(1)
                plt.plot(plot_points, dvals, color=col[i], label=desc[i])
            else:
                plt.figure(0)
                plt.plot(plot_points, vals, color=col[i])
                plt.figure(1)
                plt.plot(plot_points, dvals, color=col[i])
    plt.figure(0)
    plt.legend()
    plt.tight_layout()
    plt.figure(1)
    plt.legend()
    plt.tight_layout()
    plt.show()

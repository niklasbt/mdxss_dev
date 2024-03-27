'''Helper program to compute X-ray scattering functions

'''

# imports ####################################################################
import numpy
import tables
import os.path
from .. import LIB

# globals ####################################################################
f0_WaasKirf = {}  # dictionary to hold coefficients for WK approximation
with tables.open_file(os.path.join(LIB, 'f0_WaasKirf.h5'), 'r') as fid:
    for atom in fid.root:
        f0_WaasKirf[str(atom).split()[0][1:]] = atom.Coeff[:]


# classes ####################################################################
class ScatteringFunctions:
    def __init__(self, composition, q_grid):
        self.composition = composition
        self.q = q_grid
        self.f_array, self.weights = f_array(self.composition, self.q)
        self.f1avg = f1avg(self.f_array, self.weights)
        self.f2avg = f2avg(self.f_array, self.weights)


# functions ##################################################################
def gaussian_model(atom):
    def G_(q):
        return(
            f0_WaasKirf[atom][0] * numpy.exp(
                -f0_WaasKirf[atom][6] * (q / (4 * numpy.pi))**2
            ) +
            f0_WaasKirf[atom][1] * numpy.exp(
                -f0_WaasKirf[atom][7] * (q / (4 * numpy.pi))**2
            ) +
            f0_WaasKirf[atom][2] * numpy.exp(
                -f0_WaasKirf[atom][8] * (q / (4 * numpy.pi))**2
            ) +
            f0_WaasKirf[atom][3] * numpy.exp(
                -f0_WaasKirf[atom][9] * (q / (4 * numpy.pi))**2
            ) +
            f0_WaasKirf[atom][4] * numpy.exp(
                -f0_WaasKirf[atom][10] * (q / (4 * numpy.pi))**2
            ) +
            f0_WaasKirf[atom][5] 
        )

    return(G_)


def form_factor_x(atom, q_grid):
    f_array = numpy.zeros(q_grid.shape)
    fofq = gaussian_model(atom)
    for i in range(len(q_grid)):
        f_array[i] = fofq(q_grid[i])
    return(f_array)


def f_array(composition, q_grid):
    f_array = []
    weights = []
    for atom in composition.keys():
        weights.append(composition[atom])
        f_array.append(form_factor_x(atom, q_grid))
    weights, f_array = numpy.array(weights), numpy.array(f_array).T
    return(f_array, weights)


def f1avg(F, w):
    return(F @ w / numpy.sum(w))


def f2avg(F, w):
    return(numpy.square(F) @ w / numpy.sum(w))


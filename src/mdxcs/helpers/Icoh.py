'''Helper program to load simulated coherent scattering data

The programs 'RDFcalc.py' and 'DebyeRDF.py' compute coherent scattering
intensities from MD ensembles and write the results to an HDF5 file. This
program defines an Icoh class and loader function to read and load these data.

Basic philosophy:

Assuming we have an NVT/NpT trajectory, then the user can choose to classify
the atoms into 'types'. Each pair of types contributes a partial I(Q) to the
coherent scattering; thus the HDF5 file contains each partial which can be
summed to the total. The Icoh object can also be queried for specific partials
using the Icoh.partial_by_type() function.
'''

# imports ####################################################################
import tables
import numpy


# classes ####################################################################
class AType:
    '''Atom type class

    Attributes
    ----------
    name : str
        Label for the atom type. First two characters should be reserved for
        the element, with a final character as a label.
    resname : str
        Label for the residue to which the type belongs.
    element : str
        Element of the atom type.
    N : int
        Number of atoms in the system belonging to this type.
    '''
    def __init__(self, name, resname, element, N):
        self.name = name.decode('UTF-8')
        self.resname = resname.decode('UTF-8')
        self.element = element.decode('UTF-8')
        self.N = N


class IcohPartial:
    '''A class to hold partial coherent scattering data

    This class holds coherent scattering data for atom types alpha and beta.

    Attributes
    ----------
    alpha : AType
        Instance of AType class defining type alpha.
    beta : AType
        Instance of AType class defining type beta.
    i_d : ndarray
        Array containing the distinct scattering.
    i_s : ndarray
        Array containing the self scattering (zero-filled for alpha != beta).
    i_c : ndarray
        Total coherent scattering, i_c = i_d + i_s
    '''
    def __init__(self, alpha, beta, i_dist, i_self):
        self.alpha = alpha
        self.beta = beta
        self.i_d = i_dist
        self.i_s = i_self
        self.i_c = i_dist + i_self


class Icoh:
    '''A class to hold coherent scattering data

    This class holds coherent scattering data for a collection of atom types.

    Attributes
    ----------
    q : ndarray
        Array of Q-values
    partials : list of IcohPartial objects
        List of partials summing to the total coherent scattering. For N_xi
        types of atoms, there are N_xi + N_xi * (N_xi - 1) / 2 partials.
    types : list of str
        List of atom types (= Atype.name).
    i_d : ndarray
        Array containing the distinct scattering.
    i_s : ndarray
        Array containing the self scattering (zero-filled for alpha != beta).
    i_c : ndarray
        Total coherent scattering, i_c = i_d + i_s

    Methods
    -------
    partials_by_type(types_a, types_b):
        Returns Icoh instance built from partials where one of alpha or beta is
        selected from each of the two provided lists.
    '''
    def __init__(self, q, partials):
        self.q = q
        self.partials = partials
        types, i_dist, i_self = [], [], []
        for partial in partials:
            if partial.alpha.name == partial.beta.name:
                types.append(partial.alpha.name)
            i_dist.append(partial.i_d)
            i_self.append(partial.i_s)
        self.types = types
        self.i_d = numpy.sum(numpy.array(i_dist), axis=0)
        self.i_s = numpy.sum(numpy.array(i_self), axis=0)
        self.i_c = self.i_d + self.i_s

    def partial_by_type(self, types_a, types_b):
        '''Get partial coherent scattering by type

        Returns Icoh instance built from partials where one of alpha or beta is
        selected from each of the two provided lists.

        Parameters
        ----------
        types_a, types_b : list of str
            Provide two lists of atom types. Partials are selected according to
            the set of pairs {a, b} where a belongs to types_a and b belongs to
            types_b.

        Returns
        -------
        Icoh object constructed from the selected partials.
        '''
        # check that the proper inputs ar provided
        if not isinstance(types_a, list) or not isinstance(types_b, list):
            raise ValueError('Types must be provided as lists.')
        subpartials, target_pairs = [], []
        for type_a in types_a:
            for type_b in types_b:
                if not any(
                    pair == set((type_a, type_b)) for pair in target_pairs
                ):
                    target_pairs.append(set((type_a, type_b)))
        for partial in self.partials:
            partial_pair = set((partial.alpha.name, partial.beta.name))
            for pair in target_pairs:
                if partial_pair == pair:
                    subpartials.append(partial)
        return(Icoh(self.q, subpartials))


# functions ##################################################################
def load_icoh(fname):
    '''Load coherent scattering data from HDF5 file

    Returns an Icoh object given an HDF5 file produced by DebyeRDF.py.

    Parameters
    ----------
    fname : str
        Name of the HDF5 to load.

    Returns
    -------
    icoh : Icoh
        Icoh object.
    '''
    partials = []
    with tables.open_file(fname, 'r') as fid:
        data = fid.root.i_coh
        for datum in data:
            partials.append(
                IcohPartial(
                    AType(*datum.types[0]),
                    AType(*datum.types[1]),
                    datum.i_dist[:],
                    datum.i_self[:]
                )
            )
        icoh = Icoh(fid.root.system.q[:], partials)
    return(icoh)


def avg_icoh(icoh_list):
    '''Function to average Icoh objects'''
    # normalization for averaging
    M = len(icoh_list)
    # create base partials object
    partials_avg = icoh_list[0].partials
    for partial in partials_avg:
        partial.i_s /= M
        partial.i_d /= M
        partial.i_c /= M
    N =len(partials_avg)
    # sum all partials
    for icoh in icoh_list[1:]:
        for k in range(N):
            if \
                icoh.partials[k].alpha.name == partials_avg[k].alpha.name \
                and icoh.partials[k].beta.name == partials_avg[k].beta.name:
                partials_avg[k].i_s += icoh.partials[k].i_s / M
                partials_avg[k].i_d += icoh.partials[k].i_d / M
                partials_avg[k].i_c += icoh.partials[k].i_c / M
    return(Icoh(icoh_list[0].q, partials_avg))


def sum_icoh(icoh_list):
    '''Function to sum Icoh objects'''
    partials_sum = icoh_list[0].partials
    N = len(partials_sum)
    for icoh in icoh_list[1:]:
        for k in range(N):
            if \
                icoh.partials[k].alpha.name == partials_sum[k].alpha.name \
                and icoh.partials[k].beta.name == partials_sum[k].beta.name:
                partials_sum[k].i_s += icoh.partials[k].i_s
                partials_sum[k].i_d += icoh.partials[k].i_d
                partials_sum[k].i_c += icoh.partials[k].i_c
    return(Icoh(icoh_list[0].q, partials_sum))


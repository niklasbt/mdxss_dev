#!/usr/bin/env python
'''DebyeRDF.py: Coherent scattering from RDFs

Usage: python DebyeRDF.py <RDFs.h5> <qmin> <qmax> <qstep> <output.h5>
'''
# imports ####################################################################
import numpy as np
import tables as tb
import periodictable as pt
from scipy.integrate import trapezoid


# classes ####################################################################
class AType:
    '''Atom type class

    Holds atom type information. Only for parsing data read from HDF5 files
    produced by RDFcalc.py.
    '''
    def __init__(self, name, resname, element, N):
        self.name = name.decode('UTF-8')
        self.resname = resname.decode('UTF-8')
        self.element = element.decode('UTF-8')
        self.N = N


class RDF:
    '''Class to hold RDF data

    This class holds RDF data parsed from an HDF5 file produced by RDFcalc.py.
    '''
    def __init__(self, alpha, beta, g0, r, g, volume):
        # atom types
        self.alpha = alpha
        self.beta = beta
        # limiting behavior
        self.g0 = g0
        # cell volumes
        self.V = volume
        # g(r)
        self.r = r
        self.g = g


# functions ##################################################################
def load_rdfs(fname):
    rdf_list = []
    with tb.open_file(fname, 'r') as fid:
        rdfs = fid.root.rdfs
        for rdf in rdfs:
            rdf_list.append(RDF(
                AType(*rdf.types[0]),
                AType(*rdf.types[1]),
                rdf.g0[0],
                fid.root.system.r[:],
                rdf.g[:],
                fid.root.system.V[:]
            ))
    return(rdf_list)


def integrand(q_val, rdf):
    '''Inputs a Q-value and RDF object, computes integrand.

        Parameters
        ----------
        q_val : float
            Q-value for computation.
        rdf : rdf object
            RDF used for evalauting the integrand.

        Returns
        -------
        ndarray of size (n, )
            n is the length of the r-grid. For each frame, the integrand is,

                    (1/V) * 4 * pi * r^2 * [g_ab(r) - g0] * sin(Qr)/Qr

            This is computed for each frame, and then averaged over the
            ensemble (handling both NVT and NpT calculations).
    '''
    length = rdf.V.shape[1]
    G = (1 / rdf.V) @ (4 * np.pi * rdf.r**2 * (rdf.g - rdf.g0)) / length
    if q_val == 0:
        return(G)
    else:
        return(G * np.sin(q_val*rdf.r)/(q_val*rdf.r))


def evaluate_integral(q_grid, rdf):
    '''Compute the integral over a Q-grid.

        Parameters
        ----------
        q_grid : ndarray
            Grid of Q-values to compute the integral over.
        rdf : rdf object
            RDF used for evaluating the integral.

        Returns
        -------
        ndarray of size (n, )
            n is the size of q_grid. The function evaluates,

                i(Q) = int_{r = 0}^R{ integrand(Q) }

            where integrand is given by the `integrand` function.
    '''
    # accumulate integral values
    i = []
    for q_val in q_grid:
        func = integrand(q_val, rdf)
        i_val = trapezoid(func, rdf.r)
        i.append(i_val)
    return(np.array(i))


def compute_partial(q_grid, rdf):
    '''Compute partial coherent term for given RDF.

    Evaluates the equations,

        (1) I_self(Q) = N_a * f_a(Q)^2
        (2) I_distinct(Q) = N_a(N_b - d_ab)/V * f_a(Q) * f_b(Q) * i(Q)

    For a given pair of atom types a, b. (1) only contributes when a = b, and
    i(Q) is computed using the RDF and the function `evaluate_integral`.

        Parameters
        ----------
        q_grid : ndarray
            Grid of Q-values to compute the scattering over.
        rdf : rdf object
            RDF used for evaluating the integral.

        Returns
        -------
        id_ab: ndarray of size (n, )
            n is the size of the Q-grid. This term contains the distinct
            scattering contribution for atom types a and b.

        is_ab: ndarray of size (n, )
            n is the size of the Q-grid. This term contains the self
            scattering contribution for atom types a = b.
    '''
    # define system variables
    N_alpha, N_beta = rdf.alpha.N, rdf.beta.N
    # get form factors interpolated onto q-grid
    f_alpha = pt.elements.symbol(rdf.alpha.element).xray.f0(q_grid)
    f_beta = pt.elements.symbol(rdf.beta.element).xray.f0(q_grid)
    # compute integral
    integral = evaluate_integral(q_grid, rdf).flatten()
    # two cases: alpha = beta and alpha !=beta
    if rdf.alpha.name == rdf.beta.name:
        print(f'\t\tAtom types: {rdf.alpha.name} = {rdf.beta.name}')
        print(f'\t\tN_a = {N_alpha}')
        print(f'\t\tg_0 = {rdf.g0}')
        norm = N_alpha*(N_alpha - 1)
        id_ab = norm*f_alpha**2*integral
        print('\t...computed distinct scattering contribution')
        # in this case there is a self-scattering contribution
        is_ab = N_alpha*f_alpha**2
        print('\t...computed self-scattering contribution')
    else:
        print(f'\t\tAtom types: {rdf.alpha.name} != {rdf.beta.name}')
        print(f'\t\tN_a = {N_alpha}, N_b = {N_beta}')
        print(f'\t\tg_0 = {rdf.g0}')
        norm = N_alpha*N_beta
        # in this case we multiply by 2 to do proper counting
        id_ab = 2*norm*f_alpha*f_beta*integral
        print('\t...computed distinct scattering contribution')
        is_ab = np.zeros(len(q_grid))
    return(id_ab, is_ab)


def compute_icoh(q_grid, rdf_list):
    '''Given a system of RDFs, compute the total coherent scattering.

    RDFs to compute the total coherent scattering intensity. This function
    computes each of the self-scattering and distinct scattering terms arising
    from these RDFs.

        Parameters
        ----------
        q_grid : ndarray
            Grid of Q-values to compute the scattering over.
        rdf_list : list of rdf objects
            RDFs used for the computation.

        Returns
        -------
        id_list : list of ndarrays of size (n, )
            n is the size of the Q-grid. This term contains the distinct
            scattering contributions for all pairs of atom types a, b.

        is_list : list of ndarrays of size (n, )
            n is the size of the Q-grid. This term contains the self
            scattering contribution for atom types a, b.
    '''
    N = len(rdf_list)
    id_list, is_list = [], []
    print('\nStarting elastic scattering calculation...')
    for k in range(N):
        print(f'\n\tComputing term {k+1}/{N}...')
        id_temp, is_temp = compute_partial(
            q_grid, rdf_list[k]
        )
        id_list.append(id_temp)
        is_list.append(is_temp)
    return(id_list, is_list)


def write_output(fname, rdf_list, q_grid, id_list, is_list):
    '''Write RDFs to HDF5 file

    Parameters
    ----------
    fname : str
        Filename for output file. Format 'path/to/file/file.h5'. Obtained from
        the input JSON file.
    rdf_list : list of RDF objects
        List containing the RDF objects to be saved to file.

    Returns
    -------
    None.
    '''
    class AtomMeta(tb.IsDescription):
        '''PyTables class reference'''
        name = tb.StringCol(3, pos=0)
        resname = tb.StringCol(4, pos=1)
        element = tb.StringCol(2, pos=2)
        N = tb.Int32Col(pos=3)

    with tb.open_file(fname, 'w') as fid:
        # create group for system-wide data
        system = fid.create_group(fid.root, 'system')
        fid.create_array(system, 'q', q_grid)
        # populate rdfs
        fid.create_group(fid.root, 'i_coh')
        for k in range(len(rdf_list)):
            # one group per rdf object
            icoh_group = fid.create_group(
                fid.root.i_coh,
                f'{rdf_list[k].alpha.name}_{rdf_list[k].beta.name}'
            )
            # sub-group containing atom types
            types = fid.create_table(icoh_group, 'types', AtomMeta)
            types.row['name'] = rdf_list[k].alpha.name
            types.row['resname'] = rdf_list[k].alpha.resname
            types.row['element'] = rdf_list[k].alpha.element
            types.row['N'] = rdf_list[k].alpha.N
            types.row.append()
            types.row['name'] = rdf_list[k].beta.name
            types.row['resname'] = rdf_list[k].beta.resname
            types.row['element'] = rdf_list[k].beta.element
            types.row['N'] = rdf_list[k].beta.N
            types.row.append()
            types.flush()
            # append numerical data
            fid.create_array(icoh_group, 'i_dist', id_list[k])
            fid.create_array(icoh_group, 'i_self', is_list[k])
    return


def run(inp, qmin=0.5, qmax=25.005, qstep=0.01, output='icoh.h5'):
    # get rdfs
    print(f'\nLoading RDFs from {inp}...')
    rdfs = load_rdfs(inp)
    print(
        f'\tFound\t{int((-1 + np.sqrt(1 + 8*len(rdfs)))/2)} atom types,'
        f' {len(rdfs)} RDFs'
    )
    print('\tAtom types:\n\t\t\tResidue\tElement\tN')
    for rdf in rdfs:
        if rdf.alpha.name == rdf.beta.name:
            print(
                f'\t\t{rdf.alpha.name}:\t{rdf.alpha.resname}'
                f'\t{rdf.alpha.element}\t{rdf.alpha.N}'
            )
    print('\nSetting up calculation...')
    # define q-grid
    q = np.arange(qmin, qmax, qstep)
    print(f'\tQ_min = {q[0]:.2f} 1/A')
    print(f'\tQ_max = {q[-1]:.2f} 1/A')
    print(f'\tDQ = {qstep} 1/A')
    print(
        f'\tDimension per RDF: {rdfs[0].g.shape} ~ {np.prod(rdfs[0].g.shape)}'
    )
    # compute coherent scattering intensities
    i_dist, i_self = compute_icoh(q, rdfs)
    # write output file
    print(f"\nWriting {output}...")
    write_output(output, rdfs, q, i_dist, i_self)
    print('\n...done!')
    

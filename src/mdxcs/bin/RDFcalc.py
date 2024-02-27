#!/usr/bin/env python
'''RDFCalc.py: Compute RDFs from MD trajectories

Usage: python RDFCalc.py <input.json>
'''
# imports ####################################################################
import json
import itertools
import mdtraj as md
import numpy as np
import tables as tb


# classes ####################################################################
class AType:
    '''Atom type class

    AType class holds atom type information for RDF calculations.
    '''
    def __init__(self, name, resname, element, N, replicates):
        self.name = name
        self.resname = resname
        self.element = element
        self.N = N
        self.replicates = replicates


class RDF:
    '''RDF class

    This class initializes and holds metadata for RDF calculations. Input
    ATypes a = alpha and b = beta, the limiting behavior for g_ab(r) --> infty
    (g0 = 0 or 1), the cell volume for each frame, a list of indices for the
    ATypes (pairs), configuration for mdtraj.compute_rdf(), and the length of
    the trajectory.
    '''
    def __init__(self, alpha, beta, g0, volume, pairs, config, length):
        # config for rdf calculation
        self.config = config
        self.pairs = pairs
        array_size = np.arange(
            config['r_range'][0],
            config['r_range'][1] + config['bin_width'] / 3,
            config['bin_width']
        ).shape[0] - 1
        # atom types
        self.alpha = alpha
        self.beta = beta
        # limiting behavior
        self.g0 = g0
        # cell volumes
        self.V = volume.reshape(1, length)
        # initialize arrays
        self.r = np.zeros(array_size)
        self.g = np.zeros([length, array_size])


# functions ##################################################################
def compute_rdfs(trajectory, system, config, length):
    '''Calculate RDFs, frame-wise

    Parameters
    ----------
    trajectory : mdtraj.Trajectory
        MD trajectory.
    system : list of Atypes
        List of atom types.
    config : dict of kwargs
        Configuration for mdtraj.compute_rdf().
    length : int or 'all'
        Number of frames to include in the calculation. If 'all', the entire
        trajectory is used.

    Returns
    -------
    rdfs : list of RDF objects
        A list of RDF objects, one for each unordered pair of atom types 
        (i.e., only one of g_ab(r) = g_ba(r) is calculated. The RDF holds r, 
        g(r) where g.shape = (length, len(r)). The length of r is determined 
        by config (see RDF class).
    '''
    def prepare_rdfs():
        '''Prepare RDFs for calculation.'''
        # size of trajectory
        nonlocal length
        # volumes
        volume = trajectory[:length].unitcell_volumes * 1e3  # cubic ang.
        # get atom types and selection
        types = [atom for atom in system]
        atom_pairs = \
            [(i, i) for i in types] + \
            [x for x in itertools.combinations(types, 2)]
        rdfs = []
        for pair in atom_pairs:
            selection = trajectory.top.select_pairs(
                f"element {pair[0].element} and resname '{pair[0].resname}'",
                f"element {pair[1].element} and resname '{pair[1].resname}'"
            )
            # case 1: no replicate atoms in pair, g0 = 0
            if not pair[0].replicates and not pair[1].replicates:
                rdfs.append(RDF(
                    pair[0],
                    pair[1],
                    0,
                    volume,
                    selection,
                    config,
                    length
                ))
            # case 2: at least one pair has replicate atoms, g1 = 1
            else:
                rdfs.append(RDF(
                    pair[0],
                    pair[1],
                    1,
                    volume,
                    selection,
                    config,
                    length
                ))
        return(rdfs)

    # size of trajectory
    if length == 'all':
        length = len(trajectory)
    # initialize rdfs
    print('\tInitializing...')
    rdfs = prepare_rdfs()
    # compute rdfs
    print('\tComputing RDFs...')
    for k in range(length):
        for j in range(len(rdfs)):
            print(
                f'\t\tFrame {k+1}/{length}\tRDF: '
                f'{(j+1) + k*len(rdfs)}/{len(rdfs) * length}',
                end='\r'
            )
            if rdfs[j].pairs.shape[0] == 0:
                pass
            else:
                rdfs[j].r, rdfs[j].g[k, :] = md.compute_rdf(
                    trajectory[k],
                    rdfs[j].pairs,
                    **rdfs[j].config
                    )
            # nm -> angstroms
            rdfs[j].r *= 10
    return(rdfs)


def write_rdfs(fname, rdf_list):
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
        # find the first nontrivial r-grid (in case N = 1 for some atype)
        for rdf in rdf_list:
            if not np.array_equal(rdf.r, np.zeros(len(rdf.r))):
                fid.create_array(system, 'r', rdf.r)
                break
        fid.create_array(system, 'V', rdf_list[0].V)
        # populate rdfs
        fid.create_group(fid.root, 'rdfs')
        for rdf in rdf_list:
            # one group per rdf object
            rdf_group = fid.create_group(
                fid.root.rdfs,
                f'{rdf.alpha.name}_{rdf.beta.name}'
            )
            fid.create_array(rdf_group, 'g', rdf.g)
            fid.create_array(rdf_group, 'g0', [rdf.g0])
            # sub-group containing atom types
            types = fid.create_table(rdf_group, 'types', AtomMeta)
            types.row['name'] = rdf.alpha.name
            types.row['resname'] = rdf.alpha.resname
            types.row['element'] = rdf.alpha.element
            types.row['N'] = rdf.alpha.N
            types.row.append()
            types.row['name'] = rdf.beta.name
            types.row['resname'] = rdf.beta.resname
            types.row['element'] = rdf.beta.element
            types.row['N'] = rdf.beta.N
            types.row.append()
            types.flush()
    return


def run(atypes):
    with open(atypes, 'r') as fid:
	    inp = json.load(fid)
	# prepare calculation
    print(f"\nLoading {inp['file']}...")
    traj = md.load(inp['file'])
    print('\nPreparing system...')
    print('\tAtom types:\n\t\t\tResidue\tElement\tN')
    sys = [] 
    for atom in inp['atypes'].values():
        sys.append(AType(
            atom['name'],
            atom['resname'],
            atom['element'],
            atom['N'],
            atom['replicates']
        ))
    for atom in sys:
        print(f'\t\t{atom.name}:\t{atom.resname}\t{atom.element}\t{atom.N}')
    N_frames = inp['frames']
    if N_frames == 'all':
        print(f'\tFrames:\t{len(traj)} ({N_frames})')
    else:
        print(f'\tFrames:\t{N_frames}')
	# configuration for calculation (in nm)
    config = {
        'r_range': np.array(inp['r_range']) / 10,  # convert ang. > nm
        'bin_width': inp['bin_width'] / 10  # ibid.
    }
    print(f"\tr_max:\t{config['r_range'][1] * 10:.2f} A")
    print(f"\tdr:\t{config['bin_width'] * 10} A")
    print('\nPerforming calculation...')
    # do the calculation
    rdfs = compute_rdfs(traj, sys, config, length=N_frames)
    # write output file
    print(f"\nWriting {inp['output']}...")
    write_rdfs(inp['output'], rdfs)
    print('\n...done!\n')


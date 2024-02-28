#!/usr/bin/env python
'''RDFCalc.py: Compute RDFs from MD trajectories ---------------------------
                                                                            |
  Author: Niklas B. Thompson                                                |
                                                                            |
----------------------------------------------------------------------------
'''
# imports ####################################################################
import json
import mdtraj as md
import numpy as np


# functions ##################################################################
def parse_residues(trajectory):
    '''Parse residue information from trajectory

    Inputs an MD trajectory and parses the residue information.

    Parameters
    ----------
    trajectory : mdtraj.Trajectory
        Trajectory object read in by mdtraj.

    Returns
    -------
    residues : list
        A list with elements [mdtraj.core.topology.Residue, N] giving the
        unique residue types in the trajectory and the number of replicates.
    '''
    residues = []
    for residue in trajectory.top.residues:
        if not any(residue.name == r[0].name for r in residues):
            residues.append([residue, 1])
        else:
            for k in range(len(residues)):
                if residue.name == residues[k][0].name:
                    residues[k][1] += 1
    return(residues)


def atypes_from_residue(residue):
    '''Atom types per residue

    Returns dictionary of atom types from a given residue.

    Parameters
    ----------
    residue : mdtraj.core.topology.Residue
        Residue to parse.

    Returns
    -------
    atypes : dict
        Dictionary of atom types. The naming convention is 'symbol-resname'.
        N.B. the names should be updated manually to follow the convention
        'Symbol<x>' where 'x' is a single character identifier.
    '''
    atypes = {}
    for atom in residue[0].atoms:
        if atom.element.symbol == 'VS':
            pass
        # add new type if not present in list
        else:
            if not any(
                f'{atom.element.symbol}-{residue[0].name}'
                == key for key in atypes.keys()
            ):
                atypes[f'{atom.element.symbol}_{residue[0].name}'] = {
                    'name': f'{atom.element.symbol}_{residue[0].name}',
                    'resname': f'{residue[0].name}',
                    'element': f'{atom.element.symbol}',
                    'N': residue[1],
                    'replicates': bool(residue[1] != 1)
                }
            # else update the number of atoms
            else:
                atypes[f'{atom.element.symbol}_{residue[0].name}']['N'] \
                    += residue[1]
    return(atypes)


def parse_trajectory(fname, **kwargs):
    '''Get atom types from a trajectory

    Parameters
    ----------
    fname : MD trajectory, HDF5 format
        Filename of the trajectory to parse, loaded using mdtraj. Must be in
 	the mdtraj HDF5 format. Several standard MD trajectory formats can be 	      converted to HDF5 using the mdconvert subroutine packaged with mdtraj.

    Returns
    -------
    atypes : dict
        Dictionary of atom types. Each element in a given residue type is
        assigned as a unique type.
    '''
    print('\n##################################################################')
    print('##                                                              ##')
    print('#                  ***ENTERING ParseTraj.py***                   #')
    trajectory = md.load(fname)
    print(f'\n\tTrajectory loaded...\t{fname}')
    residues = parse_residues(trajectory)
    print('\tResidues found...')
    for residue in residues:
        print(f'\t\t{residue[0].name}')
    atypes = {}
    for residue in residues:
        atypes = {**atypes, **atypes_from_residue(residue)}
    for key in kwargs.keys():
        if key == 'name' and kwargs['name'].lower() == 'y':
            atypes = rename_types(atypes)
    if not 'r_range' in kwargs.keys():
        r_max = np.floor(np.min(trajectory.unitcell_lengths) / 2) * 10
        kwargs['r_range'] = [0, r_max] # default choice, ang.
    if not 'bin_width' in kwargs.keys():
        kwargs['bin_width'] = 0.001 # default choice, ang.
    if not 'frames' in kwargs.keys():
        kwargs['frames'] = 'all'
    if not 'output' in kwargs.keys():
        kwargs['output'] = 'rdfs.h5'
    with open('config.json', 'w') as fid:
        json.dump(
        {
                'file': fname,
                'output': kwargs['output'],
                'frames': kwargs['frames'],
                'r_range': kwargs['r_range'],
                'bin_width': kwargs['bin_width'],
                'atypes': atypes
            },
            fid,
            indent=4
        )
        print('\tFile config.json written...')
    print('\n#            ***NORMAL TERMINATION OF ParseTraj.py***            #')
    print('##                                                              ##')
    print('##################################################################\n')
    return(atypes)


def rename_types(atypes):
    '''Optional function to provide user names to atom types'''
    atypes_renamed = {}
    print('\nRequested custom names...')
    for atom in atypes:
        name = input(
            f"\tRename element {atypes[atom]['element']} in residue "
            f"{atypes[atom]['resname']}? [{atom}] "
        )
        if not name:
            atypes_renamed[atom] = atypes[atom].copy()
        else:
            atypes_renamed[name] = atypes[atom].copy()
            atypes_renamed[name]['name'] = name
    return(atypes_renamed)


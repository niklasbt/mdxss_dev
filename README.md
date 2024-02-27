# Molecular Dynamics/X-ray Coherent Scattering

This repository contains code for the Python package `mdxcs`, an acronym for `m`olecular `d`ynamics/`x`-ray `c`oherent `s`cattering. The package provides code to compute coherent X-ray scattering intensities from molecular dynamics trajectories. Please note this package is under active development.

## Theory

We consider atomistic molecular dynamics trajectories (*i.e.*, non-coarse-grained), which are typically treated as quasi-periodic systems under the minimum image convention. Most simply, such trajectories consist of a series of snapshots ('frames'), each of which is defined as a periodic cell (the 'unit cell', if you like) filled with $N$ atoms (say). We should like to compute coherent X-ray scattering intensities using the structural information encoded by such trajectories.

Assuming we are considering an isotropic scattering (e.g., a macroscopically-homogeneous system, such as a liquid), then the appropriate (atomistic) formalism is the Debye scattering equation (DSE),

```math
  I(Q) = \sum_{i=1}^N{\sum_{j=1}^N{f_j*(Q)f_i(Q)\frac{\sin{Qr_{ij}}}{Qr_{ij}}}}
``` 
where $I(Q)$ is the scattering intensity as a function of the momentum transfer, $Q$, $f_i(Q)$ is the X-ray scattering form factor for atom $i$, the double-summation is over all pairs of atoms $(i, j)$, and $r_{ij}$ is the distance between atoms (scatterers) $i$ and $j$. 

## Implementation

Under construction...

## Installation

### Dependencies

As currently built, `mdxcs` depends on:
  - itertools (*likely to change*)
  - json
  - mdtraj
  - numpy
  - periodictable
  - pytables
  - scipy


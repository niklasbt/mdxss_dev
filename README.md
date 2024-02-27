# Molecular Dynamics/X-ray Coherent Scattering

This repository contains code for the Python package `mdxcs`, an acronym for `m`olecular `d`ynamics/`x`-ray `c`oherent `s`cattering. The package provides code to compute coherent X-ray scattering intensities from molecular dynamics trajectories. Please note this package is under active development.

## Theory

We consider atomistic molecular dynamics trajectories (*i.e.*, non-coarse-grained), which are typically treated as quasi-periodic systems under the minimum image convention. Most simply, such trajectories consist of a series of snapshots ('frames'), each of which is defined as a periodic cell (the 'unit cell', if you like) filled with $N$ atoms (say). We should like to compute coherent X-ray scattering intensities using the structural information encoded by such trajectories.

Assuming we are considering an isotropic scatterer (*e.g.*, a macroscopically-homogeneous system, such as a liquid), then the appropriate (atomistic) formalism is the Debye scattering equation (DSE),

```math
  I(Q) = \sum_{i=1}^N{\sum_{j=1}^N{f_j^*(Q)f_i(Q)\frac{\sin{Qr_{ij}}}{Qr_{ij}}}}
``` 
where $I(Q)$ is the scattering intensity as a function of the momentum transfer, $Q$ (which is related to the elastic scattering angle, $2\theta$, by $Q = |\mathbf{Q}| = 4\pi\sin{\theta}/\lambda$ for X-rays of wavelength $\lambda$), $f_i(Q)$ is the X-ray scattering form factor for atom $i$, the double-summation is over all pairs of atoms $(i, j)$, and $r_{ij}$ is the distance between atoms (scatterers) $i$ and $j$. 

For isotropic systems, this equation is the most general formula to calculate coherent scattering intensities, under the independent atom approximation. As such, this code implements the DSE directly, making no approximations. The accuracy of the results is limited by (a) the size of the simulation cell (a source of finite-size errors), and (b) user-controlled parameters defining the quantization of the calculation.

Nevertheless, the code does not actually compute the DSE in the form given above. Rather, taking advangate of efficient codes to compute real-space radial distribution functions from MD simulation data, rather, we calculate,

```math
I_{\mathrm{coh},f_k}(Q) = \sum_\alpha{N_\alpha f_\alpha(Q)^2} + \\
	 \sum_\alpha{\sum_\beta{f_\alpha(Q)f_\beta(Q)\frac{N_\alpha\left(N_\beta - \delta_{\alpha\beta}\right)}{V_{f_k}}\int_0^\infty{4\pi r^2 \left(g_{\alpha\beta, f_k}(r) - g_{0,\alpha\beta} \right)\frac{\sin{Qr}}{Qr}\mathrm{d}r} }}
```

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


# imports
from mdxss.bin import ParseTraj, RDFcalc, DebyeRDF
from mdxss.helpers import Icoh
import numpy as np

# parse test trajectory
ParseTraj.parse_trajectory('water.h5', name='y')

# compute RDFS
RDFcalc.run('config.json')

# compute scattering
DebyeRDF.run('rdfs.h5')

# load scattering data
I = Icoh.load_icoh('icoh.h5')
print('\n# Calculated coherent scattering intensity:')
print('## Q (1/A), I(Q) (e.u.)')
print(np.c_[I.q, I.i_c])
print('\n')

# termination message
print('***NORMAL TERMINATION OF MDXSS***')

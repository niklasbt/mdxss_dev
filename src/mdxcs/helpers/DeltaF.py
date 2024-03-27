'''Helper program to compute differential scattering

The program takes two Icoh.Icoh objects and computes the differential
according to a user provided list of solute atom types.
'''

# imports ####################################################################
import numpy
from mdxcs.xcs.xsf import ScatteringFunctions as SF


# classes ####################################################################
class DeltaF:
    '''Differential scattering class
    
    The difference is computed according to 

        DF = (Id_A - Id_B) * q * norm

    where norm is computed from the composition of the solute,

        norm = [ N_u * < >^2_U ]**-1

    Attributes
    ----------
    XX
    '''
    def __init__(self, i_a, i_b, solute_types):
        self.q = i_a.q
        self.u_types = solute_types
        composition, solvent_types = {}, []
        for partial in i_a.partials:
            if any(partial.alpha.name == atype for atype in solute_types):
                composition[f'{partial.alpha.element}'] = partial.alpha.N
            else:
                solvent_types.append(partial.alpha.name)
            if any(partial.beta.name == atype for atype in solute_types):
                composition[f'{partial.beta.element}'] = partial.beta.N
            else:
                solvent_types.append(partial.beta.name)
        self.v_types = list(set(solvent_types))
        xsf = SF(composition, i_a.q)
        N = sum([val for val in composition.values()])
        self.norm = 1/(N * xsf.f1avg**2)
        self.DF = (i_a.i_d - i_b.i_d) * self.q * self.norm
        I_uu = i_a.partial_by_type(self.u_types, self.u_types)
        I_uv = i_a.partial_by_type(self.u_types, self.v_types)
        I_vv = i_a.partial_by_type(self.v_types, self.v_types)
        self.F_uu = I_uu.i_d * self.q * self.norm
        self.F_uv = I_uv.i_d * self.q * self.norm
        self.F_vv = I_vv.i_d * self.q * self.norm
        self.F_B = i_b.i_d * self.q * self.norm
        self.DF_vv = self.F_vv - self.F_B
    

    def excluded_solvent(self, Icoh_es):
        self.F_es = Icoh_es.i_d * self.q * self.norm
        I_ee = Icoh_es.partial_by_type(['Hd', 'Od'], ['Hd', 'Od'])
        I_evstar = Icoh_es.partial_by_type(['Hd', 'Od'], ['Hb', 'Ob'])
        self.F_ee = I_ee.i_d * self.q * self.norm
        self.F_evstar = I_evstar.i_d * self.q * self.norm
        self.DF_vvstar = self.DF_vv + self.F_es


    # TO DO: implement sine-Fourier transform
    def transform(self, qmin=1, qmax=50):
        self.r, self.DG = ft_sine(self.q, self.DF, qmin, qmax)
        tmp, self.G_uu = ft_sine(self.q, self.F_uu, qmin, qmax)
        tmp, self.G_uv = ft_sine(self.q, self.F_uv, qmin, qmax)
        tmp, self.G_B = ft_sine(self.q, self.F_B, qmin, qmax)
        tmp, self.DG_vv = ft_sine(self.q, self.DF_vv, qmin, qmax)
        if hasattr(self, 'F_es'):
            tmp, self.G_es = ft_sine(self.q, self.F_es, qmin, qmax)
            tmp, self.G_ee = ft_sine(self.q, self.F_ee, qmin, qmax)
            tmp, self.G_evstar = ft_sine(self.q, self.F_evstar, qmin, qmax)
            tmp, self.DG_vvstar = ft_sine(self.q, self.DF_vvstar, qmin, qmax)



# functions ##################################################################
# TO DO: implement sine Fourier transform (`fftftog')
def ft_sine(q, F, qmin, qmax):
    '''Compute the Fourier sine transform'''
    # find qmin index
    for k in range(len(q)):
        if q[k] > qmin:
            idxl = k
            break
    # find qmax index
    for k in range(len(q)):
        if q[k] > qmax:
            idxu = k
            break
        else:
            idxu = len(q) + 1
    # compute FT
    g, rstep = fftftog(
        numpy.concatenate([F[idxl:idxu], numpy.zeros(20480)]), 
        q[1] - q[0], 
        qmin = q[idxl]
    )
    return(numpy.arange(len(g)) * rstep, g)

# TO DO: 
def fftftog(f, qstep, qmin):
    return




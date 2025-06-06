import numpy as np
from pyscf import gto
import tomli

def check_uno(noon, thresh=1.98):
    ndb = np.count_nonzero(noon > thresh)
    nex = np.count_nonzero(noon < (2.0-thresh))
    nacto = len(noon) - ndb - nex
    return nacto, ndb, nex

chemcore_atm = [
    0,                                                                  0,
    0,  0,                                          1,  1,  1,  1,  1,  1,
    1,  1,                                          5,  5,  5,  5,  5,  5,
    5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  9,  9,  9,  9,  9,  9,
    9,  9, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 18, 18, 18, 18, 18, 18, 
   18, 18, 
           18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 23, # lanthanides
               23, 23, 23, 23, 23, 23, 23, 23, 23, 34, 34, 34, 34, 34, 34, 
   34, 34, 
           34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, # actinides
               50, 50, 50, 50, 50, 50, 50, 50, 50]
def chemcore(mol):
    core = 0
    for a in mol.atom_charges():
        core += chemcore_atm[a]
    return core

def bse_bas(bas, elem):
    import basis_set_exchange
    return gto.load(basis_set_exchange.api.get_basis(bas, elements=elem, fmt='nwchem'), elem)

def get_basis(bas, elem=None):
    if bas[:4]=='x2c-':
        return bse_bas(bas, elem)
    elif bas[:3]=='jun':
        pass
    else:
        return bas

def get_xc2(templ, omega=0.4, srdft=0.5):
    if templ == 'rsblyp':
        xc2 = f'RSH({omega},1.0,-{srdft})+{srdft}*ITYH , VWN5*0.19+ LYP*0.81'
    else:
        raise ValueError(f"Unsupported template: {templ}")
    return xc2

def toml_load(toml, domain):
    with open(toml, 'rb') as f:
        e = tomli.load(f)
        if domain in e:
            return e[domain]
        else:
            return None
        
def get_config(toml):
    return toml_load(toml, 'config')
import numpy as np
from pyscf import gto
import tomli, os

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
    if templ == 'suhf':
        xc2 = 'suhf'
    elif templ == 'rsblyp':
        xc2 = f'RSH({omega},1.0,-{srdft})+{srdft}*ITYH , VWN5*0.19+ LYP*0.81'
    elif templ == 'rsblyp1':
        xc2 = f'RSH({omega},1.0,-{srdft})+{srdft}*ITYH , VWN5*0.1425+ LYP*0.6075'
    else:
        raise ValueError(f"Unsupported template: {templ}")
    return xc2

file1 = ['./config.toml', '../config.toml']
def toml_load(toml, domain=None):
    if toml=='def':
        if os.path.exists(file1[0]):
            tomlfile = file1[0]
        elif os.path.exists(file1[1]):
            tomlfile = file1[1]
        else:
            raise ValueError('default config file not found')
    else:
        tomlfile = toml
    with open(tomlfile, 'rb') as f:
        e = tomli.load(f)
        if domain is None:
            return e
        if domain in e:
            return e[domain]
        else:
            return None
        
#def get_config(toml, domain='config'):
#    return toml_load(toml, domain)
get_config = toml_load

def _get_config_value(toml, domain, default, name):
    config = toml[domain]
    if config is None:
        raise ValueError(f"Domain '{domain}' not found in {toml}")
    name = name.lower()
    if name in config:
        value = config[name]
    elif default is not None:
        value = default
    elif 'default' in config:
        value = config['default']
    else:
        raise ValueError(f"Name '{name}' not found in domain '{domain}' of {toml}")
    return value

def get_spin(toml, name):
    return _get_config_value(toml, 'spin', 0, name)

def get_charge(toml, name):
    return _get_config_value(toml, 'charge', 0, name)

def get_config_basis(toml, name):
    return _get_config_value(toml, 'basis', None, name)

def need_x2c(bas):
    if 'x2c' in bas:
        return True
    return False

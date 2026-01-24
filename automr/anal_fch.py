from mokit.lib import gaussian, rwwfn
from automr.dump_mat import dump_mo_composition

def dump_mo_composition_fch(fchname):
    mol = gaussian.load_mol_from_fch(fchname)
    mo = gaussian.mo_fch2py(fchname)
    return dump_mo_composition(mol, mo, dump=False)

def get_noon_from_fch(fchname):
    nbf, nif = rwwfn.read_nbf_and_nif_from_fch(fchname)
    return rwwfn.read_eigenvalues_from_fch(fchname, nif, 'a')
    
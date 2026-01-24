from pyscf import gto, lib
from automr.dump_mat import *

lib.num_threads(2)

mol = gto.M(atom='''Cu 0 0 0; Cl 0 0 2.05''', basis='ccpvdz').build()
mf = mol.RHF().run()

dump_mo_composition(mol, mf.mo_coeff)
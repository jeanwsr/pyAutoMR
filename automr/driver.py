from pyphf import *
from mokit.lib import fchk as fchk_mo
from mokit.lib import fchk_uno

def run_surs(mf, xc2, conv=None, df='off', chk=None):
    mf2 = SUHF(mf)
    if df != 'off':
        mf2 = mf2.density_fit()
    if xc2=='suhf':
        mf2.dft = False
    else:
        mf2.dft=True
        #mf2.xc='0.25*HF+0.75*PBE,%s*PBE'%c
        mf2.xc=xc2
    #mf2.diis_on = False
    #mf2.diis_driver='plain'
    if conv=='slow':
        mf2.diis_start_cyc = 25
        mf2.max_cycle = 100
    elif conv=='earlydiis':
        mf2.diis_start_cyc = 3
        #mf2.diis_space = 15
        mf2.max_cycle = 100
    elif conv=='earlycrazy':
        mf2.diis_start_cyc = 0
        mf2.diis_space = 15
        mf2.max_cycle = 100
    elif 'plain' in conv:
        mf2.diis_start_cyc = 5
        mf2.diis_driver = 'plain'
        mf2.max_cycle = 100
        if 'hard' in conv:
            mf2.max_cycle = 140
            mf2.conv_tol = 5e-6
    elif conv=='damp':
        mf2.diis_start_cyc = 5
        mf2.diis_driver = 'rev1'
        mf2.diis_damp = 0.2
        mf2.max_cycle = 100
    elif conv=='slowdamp':
        mf2.diis_start_cyc = 35
        mf2.diis_driver = 'rev1'
        mf2.diis_damp = 0.2
        mf2.max_cycle = 130
    elif conv=='nodiis':
        mf2.diis_on = False
        mf2.max_cycle = 130
    if mf.mol.symmetry:
        #mf2.diis_start_cyc = 5
        mf2.max_cycle = 130
        #mf2.verbose=6
    #mf2.verbose=6
    #mf2.max_cycle = 5
    if 'lev' in conv:
        mf2.level_shift = 0.2
    if chk is not None:
        mf2.dumpchk = True
        mf2.output = get_chkname(chk)
    mf2.kernel()
    
    if xc2 == 'suhf':
        return mf2
    mf2 = mf2.to_hf()
    mf2.guesshf.mo_coeff = mf2.mo_reg
    mf2.noiter=True
    mf2.kernel()
    return mf2

def get_chkname(name):
    if isinstance(name, str):
        return name
    elif isinstance(name, tuple):
        return '.'.join(name)

def fchk(mf2, name, templ, suffix=None):
    fch = f'{name}_{templ}no.fch'
    if suffix is not None:
        fch = fch.replace('.fch', f'.{suffix}.fch')
    mf2.fchk(fch)

def fchk_uhf(mf, name, suffix, no=None, noon=None):
    fchname1 = f'{name}_{suffix}'
    fchk_mo(mf, fchname1+'_uhf.fch')
    if no is not None:
        fchk_uno(mf, fchname1+'_uno.fch', no, noon)

def run_pdft(mf2):
    for xc in ['tpbe','tblyp']:
        mf3 = supdft.PDFT(mf2, xc[1:], 'dd')
        mf3.kernel()
    for xc in ['tpbe','tblyp']:
        mf4 = supdft.PDFT(mf2, xc, 'pd')
        mf4.do_split = True
        #mf4.dump_adm = name+'.h5'
        mf4.no_thresh = 1e-4
        mf4.grids_level = 3
        mf4.kernel()
    for xc in ['ftblyp','tm06l','mc23']:
        mf4 = supdft.PDFT(mf2, xc, 'pd')
        #mf4.do_split = True
        #mf4.dump_adm = name+'.h5'
        mf4.no_thresh = 1e-4
        mf4.grids_level = 3
        mf4.kernel()

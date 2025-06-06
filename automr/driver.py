from pyphf import suscf, supdft

def run_surs(mf, xc2, conv=None, df='off'):
    mf2 = suscf.SUHF(mf)
    if df != 'off':
        mf2 = mf2.density_fit()
    mf2.dft=True
    #mf2.xc='0.25*HF+0.75*PBE,%s*PBE'%c
    mf2.xc=xc2
    #mf2.diis_on = False
    #mf2.diis_driver='plain'
    if conv=='slow':
        mf2.diis_start_cyc = 25
        mf2.max_cycle = 100
    #mf2.level_shift = 0.2
    mf2.kernel()
    
    mf2 = mf2.to_hf()
    mf2.guesshf.mo_coeff = mf2.mo_reg
    mf2.noiter=True
    mf2.kernel()
    return mf2

def run_pdft(mf2):
    for xc in ['tpbe','tblyp']:
        mf3 = supdft.PDFT(mf2, xc[1:], 'dd')
        mf3.kernel()
    for xc in ['tpbe','tblyp','tm06l','mc23']:
        mf4 = supdft.PDFT(mf2, xc, 'pd')
        #mf4.do_split = True
        #mf4.dump_adm = name+'.h5'
        mf4.no_thresh = 1e-4
        mf4.grids_level = 3
        mf4.kernel()

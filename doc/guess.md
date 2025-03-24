
# Mix
```
guess.mix(xyz, bas, charge, conv='tight')
```
Similar to Gaussian's `guess=mix`. This function perform a mix over RHF HOMO and LUMO.
Use `conv='loose'` for loose converging or `cycle=0` for no iteration after the mix.

## more pyscf interoperable
Unlike `guess.mix`, `guess._mix` allows input of `gto.Mole` or any SCF type.
```
mol = gto.M( ... )
guess._mix(mol)

mf = gto.M( ... ).RKS(xc='xxx')
guess._mix(mf)
```

# Fragment guess
Hard cases of spin-polarized singlet. Use `guess.from_frag(xyz, bas, frags, chgs, spins)`.

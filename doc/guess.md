
# Mix
```
from automr import guess
mf = guess.mix(xyz, bas)
# or
mf = guess.mix(xyz, bas, conv='tight')
```
Similar to Gaussian's `guess=mix`. This function perform a mix over RHF HOMO and LUMO.
By default, `conv='loose'` which means not fully converged.
Use `conv='tight'` for tight converging and stability loop check. Use `cycle=0` for no iteration after the mix.

## more pyscf interoperable
Unlike `guess.mix`, `guess._mix` allows input of `gto.Mole` or any SCF type.
```
mol = gto.M( ... )
mf = guess._mix(mol)

mf = gto.M( ... ).RKS(xc='xxx')
mf = guess._mix(mf)
```
This is beneficial when one wants to customize the mf object (use DFT, newton, x2c, electric field, etc.).

# Fragment guess
Hard cases of spin-polarized singlet. Use `mf = guess.from_frag(xyz, bas, frags, chgs, spins)`.
```
frags: List of List of atom numbers in each frag
chgs: List of net charge in each frag
spins: List of spin (2S, not 2S+1) in each frag
```
for example, input for N2 can be
```
mf = guess.from_frag(xyz, bas, [[0],[1]], [0,0], [3,-3])
```
This function also supports `guess._from_frag(mf_or_mol, frags, chgs, spins)`.

# Flip spin
Do a high spin SCF first, localize SOMOs, then flip to low spin. There's two flip style:

1. Flip by site
```
# start from a triplet. flip all SOMO in specified site.
mf = guess._flipspin(mol, highspin=2, flipstyle='site', site=[0])
```
TODO: set target low spin.
2. Flip by LMO
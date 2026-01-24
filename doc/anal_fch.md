# anal_fch usage

This module provides lightweight helpers for reading information from Gaussian
formatted checkpoint files (FCHK) and returning data in a Python-friendly form.

## dump_mo_composition_fch

Purpose: Load a formatted checkpoint file and return the molecular orbital
composition for all MOs using the project-wide composition routine.

Signature:

```
dump_mo_composition_fch(fchname)
```

Parameters:

- `fchname`: Path to the Gaussian `.fch` or `.fchk` file.

Returns:

- A data structure produced by `automr.dump_mat.dump_mo_composition`, holding
  the AO/MO composition information for the molecule. The function does not
  write output files (`dump=False`).

Example:

```python
from automr.anal_fch import dump_mo_composition_fch

composition = dump_mo_composition_fch("water.fch")
```

## get_noon_from_fch

Purpose: Read the occupation numbers (NOONs) from a formatted checkpoint file
and return the alpha occupation eigenvalues.

Signature:

```
get_noon_from_fch(fchname)
```

Parameters:

- `fchname`: Path to the Gaussian `.fch` or `.fchk` file.

Returns:

- A sequence of eigenvalues for the alpha spin natural orbitals.

Example:

```python
from automr.anal_fch import get_noon_from_fch

noon = get_noon_from_fch("water.fch")
```

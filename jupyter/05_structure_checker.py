# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.3
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# ### **Structure Checker**

# %% [markdown]
# For the concept of Structure Checker in general see [Structure Checker User's Guide](https://docs.chemaxon.com/latest/structure-checker_user-guide.html).

# %%
from chemaxon.io import import_mol
from chemaxon.structurechecker import StructureChecker

mol = import_mol('[NH3+]C1=CC(O)=CC=C1.O')
mol

# %%
for r in StructureChecker('solvent..moleculecharge').check(mol).results:
    print('Checker name:', r.checker_name)
    print('Atoms:       ', r.atom_indices)
    print('Bonds:       ', r.bond_indices)

# %%
StructureChecker('solvent..moleculecharge').check(mol).aggregated_colored_mol

# %%
StructureChecker('solvent->removeatom..moleculecharge->neutralize').fix(mol).fixed_mol

# %%
molecule_with_bond_errors = '''
  Mrv2401 07312513342D          

 10  9  0  0  0  0            999 V2000
   -7.5997    2.4004    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.0588    1.4786    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.1708    1.5754    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4092    1.3850    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7418    1.5754    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0370    1.4496    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3129    1.5754    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5984    1.1629    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8832    2.1457    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5984    1.9879    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  3  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
  5  6  1  0  0  0  0
  6  7  1  0  0  0  0
  7  8  1  0  0  0  0
  7  9  1  0  0  0  0
  7 10  1  0  0  0  0
M  END
'''

mol = import_mol(molecule_with_bond_errors)
StructureChecker('bondangle->clean').check(mol).aggregated_colored_mol

# %%
StructureChecker('bondangle->clean').fix(mol).fixed_mol

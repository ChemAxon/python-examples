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

# %%
from chemaxon.io import import_mol
from chemaxon.standardizer import Standardizer

toluol = import_mol('toluol')
toluol

# %%
Standardizer('aromatize').standardize(toluol)

# %%
mol2 = import_mol('[NH3+]C1=C([O-])C=CC([O-])=C1')
mol2

# %%
Standardizer('neutralize..aromatize').standardize(mol2)

# %%
mol3 = import_mol('[Na+].[O-]C(=O)C1=CC=CC=C1')
mol3

# %%
Standardizer('stripsalts').standardize(mol3)

# %%
viagra = import_mol('viagra')

Standardizer('clean3d').standardize(viagra)

# %%
Standardizer('clean3d..aromatize').standardize(viagra)

# %%
Standardizer('clean3d..aromatize..clean2d').standardize(viagra)

# %%
Standardizer('replaceatoms:queryatom=C:replaceatom=[Si]').standardize(viagra)

# %%
with open('mol_w_sgroups.mol', 'r') as file:
    mol_str = file.read()

mol_w_sgroups = import_mol(mol_str)
mol_w_sgroups

# %%
Standardizer('ungroupsgroups').standardize(mol_w_sgroups)

# %%
Standardizer('ungroupsgroups..aromatize').standardize(mol_w_sgroups)

# %% [markdown]
# Using XML configuration for standardizer:

# %%
from chemaxon.io import import_mol
from chemaxon.standardizer import Standardizer

configFile = 'path/to/my_config.xml'
with open(configFile) as config:
    standardizer = Standardizer(config)
    mol = import_mol('[H].[H]C1=C([H])C([H])=C([H])C([H])=C1[H]')
    standard_mol = standardizer.standardize(mol)

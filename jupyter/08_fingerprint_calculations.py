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
# ### **Fingerprint calculations**
#

# %% [markdown]
# Chemaxon has a number of functions that you can use to generate fingerprints.
#

# %%
from chemaxon.io import import_mol
from chemaxon.fingerprints import cfp, ecfp

mol = import_mol('CCO')

cfp = cfp(mol)
ecfp = ecfp(mol, 4, 1024)

# %% [markdown]
# These functions return `chemaxon.fingerprints.Fingerprint` objects. You can get the fingerprints
# in bytes or in binary string format.

# %%
ecfp.to_bytes()

# %%
ecfp.to_binary_string()

# %% [markdown]
# You can also calculate pharmacophore fingerprints of the molecules:
#

# %%
from chemaxon.fingerprints import pharmacophore_fp
pf = pharmacophore_fp(mol)
pf

# %% [markdown]
# This method returns a `FloatVectorFingerprint`, which contains a list of `float` values.
#
# You can also calculate _Tanimoto Dissimilarity_ for the fingerprints:

# %%
from chemaxon.fingerprints import tanimoto, ecfp, pharmacophore_fp, float_vector_tanimoto

mol = import_mol('aspirin')
mol2 = import_mol('aspirin')

ecfp1 = ecfp(mol, 4, 1024)
ecfp2 = ecfp(mol2, 4, 1024)

result1 = tanimoto(ecfp1, ecfp2)
print(result1)

mol = import_mol('aspirin')
mol2 = import_mol('coffein')

pf1 = pharmacophore_fp(mol)
pf2 = pharmacophore_fp(mol2)
result2 = float_vector_tanimoto(pf1, pf2)
print(result2)

# %% [markdown]
# Reaction fingerprints can also be claucluated using the `reaction_fp` function.

# %%
from chemaxon.fingerprints import reaction_fp

reaction_str = '[#7:1][H:2].[#6:4][S:3]([#7,#9,#17,#35,#53:7])(=[O:5])=[O:6]>>[#6:4][S:3]([#7:1])(=[O:6])=[O:5]'

reaction_mol = import_mol(reaction_str)
result_reaction = reaction_fp(reaction_mol)

result_reaction

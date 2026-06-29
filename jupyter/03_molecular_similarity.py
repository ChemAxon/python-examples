# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.4
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# ### **Finding the most similar structures to a specific target molecule**

# %%
import pandas

mols = pandas.read_table('nci1000.smiles', names=['SMILES', 'NCI_ID'])
mols

# %%
from chemaxon.io import import_mol
from chemaxon.fingerprints import ecfp, tanimoto

query = ecfp(import_mol('CC(=O)NC1=CC=C(O)C=C1'), 4, 1024) # ECFP of paracetamol

mols['SIMILARITY'] = mols.apply(lambda row: tanimoto(query, ecfp(import_mol(row['SMILES']), 4, 1024)), axis = 'columns')
mols

# %%
mols.hist('SIMILARITY', bins=50)

# %%
most_similars = mols.query('SIMILARITY > 0.5').sort_values(by=['SIMILARITY'], ascending=False)
most_similars

# %% [markdown]
# ### **The most similar structure to paracetamol from the target molecules**

# %%
import_mol(most_similars.iloc[0]['SMILES'])

# %% [markdown]
# ### **The second most similar structure**

# %%
import_mol(most_similars.iloc[1]['SMILES'])

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
# ## Molecule Search

# %%
from chemaxon.io import import_mol
from chemaxon.search import MoleculeSearch, SearchType
from chemaxon.standardizer import Standardizer

query = import_mol('C1=CC=CC=C1')
target = import_mol('OCc1ccccc1')
target

# %% [markdown]
# ### **Hit indices**

# %%
substructure_search = MoleculeSearch()
hit = substructure_search.find(query, target)
hit.hit_indices

# %% [markdown]
# ### **Duplicate search, setting search options**

# %% [markdown]
# **Note:** In case of duplicate search, the used tautomerization mode is **_generic_**. Which means, that all the theoretically possible tautomers of the target  is matched with the query structure itself.

# %%
duplicate_search = MoleculeSearch(search_type=SearchType.DUPLICATE)
hit = duplicate_search.find(query, target)
hit is None

# %% [markdown]
# ### **Hit coloring**

# %%
q = import_mol('CCCCCC')
t = import_mol('C1CCCCC1Cl')
MoleculeSearch().find(q, t, return_colored_hit=True).colored_hit

# %% [markdown]
# ### **Setting standardizer**

# %%
basic_arom_search = MoleculeSearch(standardizer=Standardizer('aromatize:b'))
q_basic = import_mol('CC')
t_basic = import_mol('O=C1NC=CC=C1')
t_basic

# %%
basic_arom_search.find(q_basic,t_basic, return_colored_hit=True).colored_hit

# %% [markdown]
# ### **No standardizer**

# %%
search_no_st = MoleculeSearch(standardizer="")
search_no_st.find(import_mol('c1ccccc1'), import_mol('c1ccccc1CC'), True).colored_hit

# %%
search_no_st = MoleculeSearch(standardizer=None)
search_no_st.find(import_mol('c1ccccc1'), import_mol('C1=CC=CC=C1'), True) is None

# %% [markdown]
# ### **Bulk Search**

# %%
bulk_search = MoleculeSearch()

query = import_mol('C1=CC=CC=C1')
targets = [import_mol('OCc1ccccc1'), import_mol('CC'), import_mol('c1ccccc1')]
hits = bulk_search.find_in_list(query, targets)
hits

# %%
[hit.hit_indices if hit is not None else hit for hit in hits]

# %%
from chemaxon.io import export_mol
colored_hits = bulk_search.find_in_list(query, targets, return_colored_hit=True)
colored_hits[0].colored_hit

# %% [markdown]
# ### **Bulk search speed**

# %%
import time 
start = time.time()
hits = bulk_search.find_in_list(query, targets*5000)
time.time()-start

# %%
substructure_search = MoleculeSearch()
start = time.time()
hits = [substructure_search.find(query, t) for  t in targets*5000 ]
time.time()-start

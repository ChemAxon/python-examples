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
# # __Pandas examples with Chemaxon molecules__

# %% [markdown]
# You can see simpler examples of using `pandas.DataFrame.read_table` with Chemaxon molecules in the [Calculators](02_calculators.ipynb) and [Molecular similiarity notebook](03_molecular_similarity.ipynb). But here are some additional examples showing how to work with Chemaxon molecules in [pandas.DataFrames](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html).

# %% [markdown]
# By default, molecules are being printed in the _cxsmiles_ format. In case of issues with that, it is written in the _cxsmarts_ format. If none of them are possible for some reason, representation falls back to use the original format, which was recognized during the import process.

# %%
import pandas as pd
from chemaxon.io import import_mol, open_for_import

mol = import_mol('CC(=O)NC1=CC=C(O)C=C1')

mol_lst = []
with open_for_import('mol_w_sgroups.mol') as mol_importer:
    for m in mol_importer:
        mol_lst.append(m)

d = {'molecule': [mol] + mol_lst }
df_mols = pd.DataFrame(data=d)

df_mols

# %% [markdown]
# We also provide a utility function that returns a dictionary object that already contains the `Molecule` objects, so you can easily create a `DataFrame` from it.

# %%
from chemaxon.pandasutil import load_molecules_for_pandas, prepare_molecules_for_pandas

df_mols_2 = pd.DataFrame(data=load_molecules_for_pandas('mol_w_sgroups.mol'))
df_mols_2

# %% [markdown]
# Utility method `prepare_molecules_for_pandas` enables creating a `DataFrame` object from a list of `chemaxon.Molecule` objects. `read_properties_to_columns` parameter facilitates the creation of columns from _sdf_ molecule properties.

# %%
from chemaxon.pandasutil import prepare_molecules_for_pandas

mol_lst_with_properties = []
with open_for_import('mol_with_properties.sdf') as mol_importer:
    for m in mol_importer:
        mol_lst_with_properties.append(m)

df_mols_with_properties = pd.DataFrame(data=prepare_molecules_for_pandas(mol_lst_with_properties, read_properties_to_columns=True))
df_mols_with_properties

# %% [markdown]
# In case of creating HTML output from a `pandas.DataFrame`, you can use the helper function `mol_to_svg_formatter` to visualize the molecules as _SVG_ images.

# %%
from chemaxon.pandasutil import mol_to_svg_formatter

df_mols.to_html('web_view.html', escape=False, formatters=dict(molecule=mol_to_svg_formatter))

# %% [markdown]
# HTML output combined with molecules loaded from an input file. File load is parametrized in order to see the exported molecule in the _cxsmiles_ column:

# %%
from chemaxon.pandasutil import mol_to_svg_formatter
df_mol_with_custom_columns = pd.DataFrame(data=load_molecules_for_pandas(file_path='mol_w_sgroups.mol', mol_obj_column="mols", mol_str_column="cxsmiles"))
df_mol_with_custom_columns.to_html('web_view_cxsmiles_column.html', escape=False, formatters=dict(mols=mol_to_svg_formatter))

# %% [markdown]
# Since the `Molecule` objects are being stored in the `DataFrame`, not just their representation, you can easily calculate properties for them and store the results in new columns.

# %%
from chemaxon.calculations import logp

df_mols['LogP'] = df_mols['molecule'].apply(lambda m: logp(m))
df_mols

# %% [markdown]
# You can also easily create new molecule columns based on existing columns, that contain molecules in any of the supported formats.

# %%
d = {'SMILES': ['CN1C=NC2=C1C(=O)N(C)C(=O)N2C'], 'name': ['coffein'] }
df = pd.DataFrame(data=d)

df = pd.concat([df, pd.DataFrame.from_records([{'SMILES' : 'CC[C@H](C)[C@@H]1NC(=O)[C@H](CC2=CC=C(O)C=C2)NC(=O)[C@@H](N)CSSC[C@H](NC(=O)[C@H](CC(N)=O)NC(=O)[C@H](CCC(N)=O)NC1=O)C(=O)N1CCC[C@H]1C(=O)N[C@@H](CC(C)C)C(=O)NCC(N)=O', 'name' : 'oxytocin'}])])

df['molecule'] = df['SMILES'].apply(lambda s: import_mol(s))
df


# %% [markdown]
# You can also use molecule properties in `DataFrame` objects:

# %%
from chemaxon.io import open_for_import

with open_for_import('mol_with_properties.sdf') as mol_iterator:
    mols = list(mol_iterator)

d = {'molecule': mols }
df_props = pd.DataFrame(data=d)
df_props['string property'] = df_props['molecule'].apply(lambda m: m.get_property('test_str_property'))
df_props['int array property'] = df_props['molecule'].apply(lambda m: m.get_property('test_int_array_property'))

df_props

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
# ### **Installing the Chemaxon Python module**
#
# Install guide: https://docs.chemaxon.com/display/docs/python-api_installation.md

 # %%
 import sys

 # !{sys.executable} -m pip install chemaxon

# %% [markdown]
# ### **Checking the details of the installed Python API license**

# %%
[lic for lic in chemaxon.licenses() if lic['product'] == 'Python API']

# %% [markdown]
# ### **Molecule import / export**

# %%
from chemaxon.io import import_mol, export_mol
mol = import_mol('CC(=O)NC1=CC=C(O)C=C1')
mol

# %%
print(export_mol(mol, 'smiles:u'))

# %%
print(export_mol(mol, 'smiles:+H'))

# %%
print(export_mol(mol, 'mol'))

# %% [markdown]
# ### **Reading multiple molecules from file**

# %%
from chemaxon.io import open_for_import

with open_for_import('nci1000.smiles') as mol_iterator:
    mols = list(mol_iterator)

len(mols)

# %%
mols[0]

# %%
mols[598]

# %%
from chemaxon.io import open_for_export, open_for_import

with open_for_export('nci1000.mrv', 'mrv') as exporter:
    export_res = all(exporter.write(m) for m in mols)

# test exported file:
with open_for_import('nci1000.mrv') as importer:
    mol_reimp = list(importer)

len(mol_reimp)

# %% [markdown]
# ### **Type checking of function parameters**

# %%
from chemaxon.io import import_mol

import_mol('c1ccccc1', True) # erroneous parametrization raises TypeError

# %% [markdown]
# ### **Molecule property handling**

# %% [markdown]
# You can store custom properties in the __Molecule__ objects using the __set_property__ and __get_property__ methods. A dictionary of the __Molecule__ object's properties can also be reteived using the __get_properties_dict__ method.
#
# _**Note:** Modifying that dictionary will not change the properties of the __Molecule__ object. For that, you should use __set_property__ method._

# %% [markdown]
# Supported property value types are: __str__, __int__, __float__, __bool__, __list__ (of int, float). Considering the list of scalar values, the api is aligned to the behavior of Chemaxon Java API. So, double values can only be used with dots, not with commas.

# %% [markdown]
# #### Reading properties:

# %%
from chemaxon.io import open_for_import

with open_for_import('mol_with_properties.sdf') as mol_iterator:
    mols = list(mol_iterator)

mols[0]

# %%
from chemaxon.io import open_for_import

with open_for_import('mol_with_properties.sdf') as mol_iterator:
    mols = list(mol_iterator)

mol = mols[0]

print('Number of properties: ' + str(len(mol.get_property_dict())))
print('\ttest_str_property -> ' + mol.get_property('test_str_property'))
print('\ttest_double_property -> ' + str(mol.get_property('test_double_property')))
print('\ttest_int_property -> ' + str(mol.get_property('test_int_property')))
print('\ttest_boolean_property -> ' + str(mol.get_property('test_boolean_property')))
print('\ttest_int_array_property -> ' + str(mol.get_property('test_int_array_property')))
print('\ttest_double_array_property -> ' + str(mol.get_property('test_double_array_property')))

print()
print('All properties as dict: \n' + str(mol.get_property_dict()))

# %% [markdown]
# #### Adding properties to the molecule object
#
# Supported types are: __bool__, __int__, __float__, __str__, __list__ of int, __list__ of float__. Note that the list of scalar values should be provided as a list of the corresponding type, not as a string. So, double values can only be used with dots, not with commas.

# %%
from chemaxon.io import import_mol

mol = import_mol('CC(=O)NC1=CC=C(O)C=C1')

mol.set_property('test_str_property', 'test string')
mol.set_property('test_double_property', 3.14)
mol.set_property('test_int_property', 42)
mol.set_property('test_boolean_property', True)
mol.set_property('test_int_array_property', [1, 2, 3])
mol.set_property('test_double_array_property', [1.1, 2.2, 3.3])

print('Number of properties: ' + str(len(mol.get_property_dict())))
print('\ttest_str_property -> ' + mol.get_property('test_str_property'))


# Properties are persisted on export of course:
from chemaxon.io import export_mol

print(export_mol(mol, 'sdf'))

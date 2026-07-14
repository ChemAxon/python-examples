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
# ### **Calculating physico-chemical properties**
#
# Currently available calculations: https://apidocs.chemaxon.com/python_api/apidocs/chemaxon/calculations.html

# %%
from chemaxon.io import import_mol
from chemaxon.calculations import logp, logd, pka, hlb

mol = import_mol('CC(=O)NC1=CC=C(O)C=C1') # paracetamol

print('logP:        ', logp(mol))
print('logD[pH 9.0]:', logd(mol, ph=9.0))
print('pKa:         ', pka(mol))
print('hlb:         ', hlb(mol))

# %%
from chemaxon.io import import_mol
from chemaxon.calculations import pka

mol = import_mol('aspirin')

pka_result = pka(mol)
pka_result.mol

# %%
import sys

# !{sys.executable} -m pip install matplotlib

# %%
from chemaxon.calculations import logd_ph_range, PhRange
import matplotlib.pyplot as pyplt

res = logd_ph_range(import_mol('aspirin'), PhRange(0, 14, 0.5))
pyplt.plot([r.ph for r in res],[r.logd for r in res])
pyplt.title('logD by pH')
pyplt.xlabel('pH')
pyplt.ylabel('logD')

# %% [markdown]
# ### **Chemical Terms**
#
# Please note that Chemaxon Python API **does not support** all available [Chemical Terms](https://docs.chemaxon.com/display/docs/chemical-terms_functions-by-categories.md) functions. Check the documentation for details.

# %%
from chemaxon.calculations import evaluate

print('Formula:   ', evaluate(mol, 'formula()'))
print('Atom count:', evaluate(mol, 'atomCount()'))
print('Ring count:', evaluate(mol, 'ringCount()'))
print('Fsp3:      ', evaluate(mol, 'fsp3()'))
print('logP:      ', evaluate(mol, 'logP()'))
print('logS:      ', evaluate(mol, 'logS()'))

# %% [markdown]
# #### **Tautomer region of warfarin**

# %%
warfarin = import_mol('warfarin')
import_mol(evaluate(warfarin, 'molFormat(genericTautomer(), "smarts")'))

# %% [markdown]
# ### **Calculating properties for a set of molecules**

# %%
import sys

# !{sys.executable} -m pip install pandas matplotlib

# %%
import pandas

mols = pandas.read_table('nci1000.smiles', names=['SMILES', 'NCI_ID'])
mols

# %%
mols['LOGP'] = mols.apply(lambda row: round(logp(import_mol(row['SMILES'])), 2), axis = 'columns')
mols['LOGD[3.0]'] = mols.apply(lambda row: logd(import_mol(row['SMILES']), ph=3.0), axis = 'columns')
mols['LOGD[7.4]'] = mols.apply(lambda row: logd(import_mol(row['SMILES']), ph=7.4), axis = 'columns')
mols['LOGD[11.0]'] = mols.apply(lambda row: logd(import_mol(row['SMILES']), ph=11.0), axis = 'columns')
mols

# %%
mols.hist(['LOGP', 'LOGD[3.0]', 'LOGD[7.4]', 'LOGD[11.0]'], bins=50)

# %% [markdown]
# ### **ADMET predictors**

# %% [markdown]
# There are four endpoints available in the Chemaxon Python API falling into the ADMET predictions category: _herg_classifiaction_, _herg_activity_, _bbb_ (Blood-Brain Barrier) and _cns_mpo_ (CNS Multiparameter Optimisation).

# %%
from chemaxon.calculations import herg_classification, herg_activity, bbb, cns_mpo

print(str(mol), ' (aspirin)')
print()
print('HERG classification:', 'SAFE' if herg_classification(mol).classification == 0 else 'TOXIC')
print('HERG activity:      ', herg_activity(mol).value)
print('BBB score:          ', bbb(mol).score)
print('CNS MPO score:      ', cns_mpo(mol).score)

# %%
mols_admet = pandas.read_table('nci50.smiles', names=['SMILES', 'NCI_ID'])

mols_admet['HERG classification'] = mols_admet.apply(lambda row: 'SAFE' if herg_classification(import_mol(row['SMILES'])).classification == 0 else 'TOXIC', axis = 'columns')
mols_admet['HERG activity'] = mols_admet.apply(lambda row: herg_activity(import_mol(row['SMILES'])).value, axis = 'columns')
mols_admet['BBB score'] = mols_admet.apply(lambda row: bbb(import_mol(row['SMILES'])).score, axis = 'columns')
mols_admet['CNS MPO score'] = mols_admet.apply(lambda row: cns_mpo(import_mol(row['SMILES'])).score, axis = 'columns')

mols_admet

# %% [markdown]
# Filter out the TOXIC compounds based on their hERG classification.

# %%
mols_admet.filter(items=['SMILES', 'NCI_ID', 'HERG classification', 'HERG activity']).query('`HERG classification` == "TOXIC"')

# %% [markdown]
# ### **Isoelectric point calculation**

# %%
from chemaxon.calculations import isoelectric_point

glycine_mol = import_mol('C(C(=O)O)N')

res = isoelectric_point(glycine_mol)

print('Isoelectric point of glycine: ', res.isoelectric_point)
print('Charge distributions: \n', res.charge_distributions)
print()

# %%
from chemaxon.calculations import PhRange

res = isoelectric_point(glycine_mol, ph_range=PhRange(0, 14, 0.5))

rows = []
for charge_res in res.charge_distributions:
    rows.append({'pH': charge_res.ph, 'Charge distribution': charge_res.charge})

df_isoelectric = pandas.DataFrame(rows)
df_isoelectric

# %%
# #!{sys.executable} -m pip install matplotlib

import matplotlib

df_isoelectric.plot(x='pH', y='Charge distribution', kind='line', marker='o', title='Charge Distribution vs pH for Glycine', grid=True)


# %% [markdown]
# ### **Conformer calculation**

# %%
#helper function to show the result of the conformer calculation
def display_result(result: list):
    for res in result:
        display(res[0])
        print(res[1])
        
from chemaxon.io import import_mol
from chemaxon.calculations import conformers

mol = import_mol('OC1CC(O)CCC1')
result = conformers(mol)
display_result(result)

# %%
from chemaxon.calculations import ConformerOptions, ConformerForceField, ConformerEnergyUnit, ConformerOptimization

options = ConformerOptions()
options.max_number_of_conformers = 3
options.energy_unit = ConformerEnergyUnit.KJ_PER_MOL
options.force_field = ConformerForceField.MMFF94
options.optimization_limit = ConformerOptimization.VERY_STRICT
result = conformers(mol, options=options)
display_result(result)

# %% [markdown]
# ### **Solubility calculation**

# %%
from chemaxon.io import import_mol
from chemaxon.calculations import solubility

mol_str = 'CC(=O)OC1=CC=CC=C1C(O)=O'
mol = import_mol(mol_str)
solubility(mol)

# %%
from chemaxon.calculations import solubility_ph_range, PhRange, SolubilityUnit

result = solubility_ph_range(mol, PhRange(0.0, 14.0, 1.0), unit = SolubilityUnit.MOL_PER_L)
rows = []
for res in result:
    rows.append({'pH': res.ph, 'Solubility': res.solubility})

df_solubility = pandas.DataFrame(rows)
df_solubility

# %% [markdown]
# ### **Tautomer calculation**

# %%
from chemaxon.io import import_mol 
from chemaxon.calculations import (all_tautomers, dominant_tautomer_distribution, 
    canonical_tautomer, major_tautomer, TautomerAdvancedOptions)

mol = import_mol('OC1=NC=CC2=CC=NC=C12')
result = all_tautomers(mol)
for res in result:
    display(res)


# %%
#helper function to show the result of tautomer calculation
def display_tautomer_result(result: list):
    for res in result:
        display(res[0])
        print(res[1])

mol = import_mol('CC1=CNCC(O)=C1')
result = dominant_tautomer_distribution(mol)
display_tautomer_result(result)

# %%
mol = import_mol('OC1N(S)CC=CC1=O')
canonical_tautomer(mol, normal=True)

# %%
mol = import_mol('OC1=NC=CC2=CC=NC=C12')
major_tautomer(mol, ph=14.0)

# %%
mol = import_mol('CC(=O)\\C=C(/O)CC1=CC=CC=C1')
advanced_options = TautomerAdvancedOptions()
advanced_options.protect_double_bond_stereo = True
major_tautomer(mol, options=advanced_options)

# %% [markdown]
# ### **Geometrical calculations**

# %% [markdown]
# #### Polar Surface Area (2D) 

# %% [markdown]
# For more information about this plugin, check the [public documentation](https://docs.chemaxon.com/latest/calculators_polar-surface-area-plugin-2d.html).

# %%
from chemaxon.calculations import polar_surface_area

mol_psa = import_mol('CC(=O)OC1=CC=CC=C1C(O)=O')
psa_value = polar_surface_area(mol_psa)

print('The molecule: (aspirin)')
display(mol_psa)
print('Polar Surface Area value: ' + str(psa_value))


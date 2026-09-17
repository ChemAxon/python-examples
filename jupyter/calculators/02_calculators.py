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

# %% [markdown]
# ### **Visualizing atom-index-based results**
#
# Many calculations return a value for every atom (e.g. `pka_result.pka_values` above).
# Function `visualize_atom_values` renders such per-atom results as an SVG image placing each
# value next to its atom. Required parameters are the molecule and the result list. Since `pka_values`
# isn't rounded, it can be rounded by setting the `precision` parameter (default value is `2`).

# %%
from IPython.display import SVG
from chemaxon.io import visualize_atom_values
from chemaxon.calculations import AtomDoubleValue

SVG(visualize_atom_values(mol,  pka_result.pka_values, precision=2))

# %% [markdown]
# When the result list carries more than one value per atom (e.g. `ChargeValue` has both
# `formal_charge` and `total_charge`), the `value_attribute` parameter should be used to pick the one to display:

# %%
from chemaxon.calculations import charge_by_atoms, ChargeValue

charge_result = charge_by_atoms(mol)
SVG(visualize_atom_values(mol, charge_result.charge_values, value_attribute='total_charge', precision=3))

# %% [markdown]
# Visualized result can be exported to svg file as well using standard python functions: 

# %%
svg_content = visualize_atom_values(mol, charge_result.charge_values, value_attribute='total_charge', precision=3)
with open("aspirinCharge.svg", "w") as f: f.write(svg_content)

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

mols = pandas.read_table('../nci1000.smiles', names=['SMILES', 'NCI_ID'])
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
mols_admet = pandas.read_table('../nci50.smiles', names=['SMILES', 'NCI_ID'])

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
from chemaxon.calculations import ConformerOptions, ConformerForceField, EnergyUnit, OptimizationLimit

options = ConformerOptions()
options.max_number_of_conformers = 3
options.energy_unit = EnergyUnit.KJ_PER_MOL
options.force_field = ConformerForceField.MMFF94
options.optimization_limit = OptimizationLimit.VERY_STRICT
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
# ### **Resonance calculation**
#
# The resonance plugin enumerates the resonant (mesomeric) structures of a molecule — the set of
# Lewis structures that differ only in the distribution of the electrons, not in the position of the
# atoms. By default `resonance_structures` returns only the *major contributors*, filters out
# symmetrical duplicates and generates at most 1000 structures.

# %%
from chemaxon.io import import_mol
from chemaxon.calculations import resonant_structures, canonical_resonant_structure

mol = import_mol('[O-]C(=O)C1=CC=CC=C1')
result = resonant_structures(mol)
for res in result:
    display(res)

# %% [markdown]
# Pass `major_contributors_only=False` to obtain every contributor, and `symmetry_filtering=False`
# to keep symmetrical duplicates. `max_structures` caps the number of generated structures.

# %%
mol = import_mol('[O-]C(=O)C([O-])=O')
all_contributors = resonant_structures(mol, major_contributors_only=False, symmetry_filtering=False)
major_contributors = resonant_structures(mol)
print(f'{len(all_contributors)} contributors in total, {len(major_contributors)} of them major')

# %% [markdown]
# `canonical_resonant_structure` returns a single, canonical representative form, which is useful
# for example for structure indexing and searching. Set `clean_structure=True` to lay the result
# out in 2D.

# %%
mol = import_mol('[O-]C(=O)C1=CC=CC=C1')
canonical_resonant_structure(mol, clean_structure=True)

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


# %% [markdown]
# ### **Hückel Analysis**

# %% [markdown]
# Several methods are available that are based on the Hückel Molecular Orbital (HMO) theorem/method. An example for each on a `chemaxon.Molecule` (these can also be calculated on specific pH values):

# %%
from IPython.display import SVG
from chemaxon.calculations import (hmo_electrophilic_localization_energy,
                                    hmo_nucleophilic_localization_energy, hmo_electron_density,
                                    hmo_charge_density, hmo_electrophilic_order,
                                    hmo_nucleophilic_order, hmo_pi_energy)
from chemaxon.io import import_mol, visualize_atom_values

mol_hmo = import_mol('CN1C=NC2=C1C(=O)NC(=O)N2C')

electr_localization_energy = hmo_electrophilic_localization_energy(mol_hmo)
nucleo_localization_energy = hmo_nucleophilic_localization_energy(mol_hmo)
electr_dens = hmo_electron_density(mol_hmo)
ch_dens = hmo_charge_density(mol_hmo)
electr_ord = hmo_electrophilic_order(mol_hmo)
nucl_ord = hmo_nucleophilic_order(mol_hmo)
pi_energy = hmo_pi_energy(mol_hmo)

print('Pi energy: '+ str(pi_energy))

# %%
print('L(+) localization energy by atom indices:')
SVG(visualize_atom_values(mol_hmo, electr_localization_energy))

# %%
print('L(-) localization energy by atom indices:')
SVG(visualize_atom_values(mol_hmo, nucleo_localization_energy))

# %%
print('Electrophilic densities by atom indices:')
SVG(visualize_atom_values(mol_hmo, electr_dens))

# %%
print('Charge densities by atom indices:')
SVG(visualize_atom_values(mol_hmo, ch_dens))

# %%
print('E(+) orders by atom indices:')
SVG(visualize_atom_values(mol_hmo, electr_ord))

# %%
print('Nu(-) orders by atom indices:')
SVG(visualize_atom_values(mol_hmo, nucl_ord))

# %% [markdown]
# ### **Molecular Surface Area 3D**

# %% [markdown]
# There are two versions available for 3D molecular surface area calculation:
# - [van der Waals molecular surface area](https://docs.chemaxon.com/latest/calculators_molecular-surface-area-plugin-3d.html)
# - [solvent accessible molecular surface area](https://docs.chemaxon.com/latest/calculators_molecular-surface-area-plugin-3d.html)

# %%
from chemaxon.io import import_mol, export_mol
from chemaxon.calculations import van_der_waals_surface_area, solvent_accessible_surface_area

aspirin = import_mol('CC(=O)OC1=CC=CC=C1C(O)=O')
vdw_result = van_der_waals_surface_area(aspirin)
asa_result = solvent_accessible_surface_area(aspirin)

print('van der Waals surface area of aspirin: ' + str(vdw_result.surface_area))
print('ASA surface area of aspirin: ' + str(asa_result.surface_area))
print('ASA of atoms with positive partial charge: ' + str(asa_result.asa_plus))
print('ASA of atoms with negative partial charge: ' + str(asa_result.asa_negative))
print('ASA of hydrophobic atoms (|partial charge| < 0.125): ' + str(asa_result.asa_hydrophobic))
print('ASA of polar atoms (|partial charge| >= 0.125): ' + str(asa_result.asa_polar))

# %%
from IPython.display import SVG
from chemaxon.io import visualize_atom_values

print('van der Waals increments of aspirin:')
SVG(visualize_atom_values(aspirin, vdw_result.increments))

# %% [markdown]
# ### **Refractivity calculation**
#
# Molar refractivity is a descriptor of the molecular volume and the London dispersive forces playing a role in drug-receptor interactions. For more information, check the [public documentation](https://docs.chemaxon.com/latest/calculators_refractivity-plugin.html).

# %%
from chemaxon.io import import_mol
from chemaxon.calculations import refractivity

mol_refr = import_mol('CC(=O)OC1=CC=CC=C1C(O)=O')  # aspirin

result = refractivity(mol_refr)

print('The molecule: (aspirin)')
display(mol_refr)
print('Molar refractivity:', result.molar_refractivity)

# %%
from IPython.display import SVG
from chemaxon.io import visualize_atom_values

print('Atomic refractivity increments:')
SVG(visualize_atom_values(mol_refr, result.refractivity_values, value_attribute='refractivity'))

# %%
from IPython.display import SVG
from chemaxon.io import visualize_atom_values

print('Atomic hydrogen refractivity increments:')
SVG(visualize_atom_values(mol_refr, result.refractivity_values, value_attribute='hydrogen_refractivity'))

# %% [markdown]
# ### **Geometrical Descriptors calculation**
#
# Molecule-scope 3D descriptors of a conformation (Dreiding and MMFF94 strain energy, minimal/maximal projection area and radius, minZ/maxZ, van der Waals volume), plus atom-level distance/angle/dihedral measurements and per-atom steric hindrance. A 3D conformer is generated automatically for input without 3D coordinates. For more information, check the [public documentation](https://docs.chemaxon.com/latest/calculators_geometrical-descriptors-plugin.html).

# %%
from chemaxon.io import import_mol
from chemaxon.calculations import geometrical_descriptors

mol_geom = import_mol('CC(=O)OC1=CC=CC=C1C(O)=O')  # aspirin

result = geometrical_descriptors(mol_geom)

print('The molecule: (aspirin)')
display(mol_geom)
print('Dreiding energy:', result.dreiding_energy)
print('MMFF94 energy:', result.mmff94_energy)
print('Minimal projection area:', result.minimal_projection_area)
print('Maximal projection area:', result.maximal_projection_area)
print('Volume:', result.volume)

# %% [markdown]
# The plugin also exposes standalone atom-level measurements (distance, angle, dihedral) and steric hindrance:

# %%
from chemaxon.calculations import distance, angle, dihedral, steric_hindrance

print('Distance between atom 0 and atom 1:', distance(mol_geom, 0, 1))
print('Angle at atom 1 (atoms 0-1-2):', angle(mol_geom, 0, 1, 2))
print('Dihedral of atoms 0-1-2-3:', dihedral(mol_geom, 0, 1, 2, 3))

# %%
from IPython.display import SVG
from chemaxon.io import visualize_atom_values

print('Steric hindrance per atom:')
SVG(visualize_atom_values(mol_geom, steric_hindrance(mol_geom)))

# %% [markdown]
# ### **Structural Frameworks calculation**
#
# Reduces a molecule to a structural framework (scaffold), such as the Bemis-Murcko scaffold or a ring system, by stripping side chains, generalizing atoms/bonds or selecting ring systems. For more information, check the [public documentation](https://docs.chemaxon.com/display/docs/calculators_structural-frameworks-plugin.html).

# %%
from chemaxon.io import import_mol
from chemaxon.calculations import structural_framework, FrameworkType

mol_scaffold = import_mol('CC(=O)OC1=CC=CC=C1C(O)=O')  # aspirin

print('The molecule: (aspirin)')
display(mol_scaffold)

bemis_murcko = structural_framework(mol_scaffold, framework_type=FrameworkType.BEMIS_MURCKO)
print('Bemis-Murcko framework:')
display(bemis_murcko)

ring_systems = structural_framework(mol_scaffold, framework_type=FrameworkType.ALL_RING_SYSTEMS)
print('All ring systems:')
display(ring_systems)

# %% [markdown]
# ### **Hydrogen Bond Donor/Acceptor (HBDA) calculation**
#
# For more information, check the [public documentation](https://docs.chemaxon.com/latest/calculators_hydrogen-bond-donor-acceptor-plugin.html).

# %%
from chemaxon.io import import_mol
from chemaxon.calculations import hbda

mol_hbda = import_mol('CC(=O)NC1=CC=C(O)C=C1')  # paracetamol

result = hbda(mol_hbda)

print('The molecule: (paracetamol)')
display(mol_hbda)
print('Donor atom count:', result.donor_atom_count)
print('Acceptor atom count:', result.acceptor_atom_count)
print('Donor site count:', result.donor_site_count)
print('Acceptor site count:', result.acceptor_site_count)

# %%
from IPython.display import SVG
from chemaxon.io import visualize_atom_values

print('Per-atom donor counts:')
SVG(visualize_atom_values(mol_hbda, result.atom_values, value_attribute='donor_count'))

# %%
from IPython.display import SVG
from chemaxon.io import visualize_atom_values

print('Per-atom acceptor counts:')
SVG(visualize_atom_values(mol_hbda, result.atom_values, value_attribute='acceptor_count'))

# %% [markdown]
# The donor/acceptor counts can also be calculated for the major microspecies over a pH range:

# %%
from chemaxon.calculations import hbda_ph_range, PhRange
import matplotlib.pyplot as pyplt

ph_results = hbda_ph_range(mol_hbda, PhRange(0, 14, 0.5))
pyplt.plot([r.ph for r in ph_results], [r.donor_count for r in ph_results], label='Donor count')
pyplt.plot([r.ph for r in ph_results], [r.acceptor_count for r in ph_results], label='Acceptor count')
pyplt.title('HBDA counts by pH (paracetamol)')
pyplt.xlabel('pH')
pyplt.ylabel('Count')
pyplt.legend()

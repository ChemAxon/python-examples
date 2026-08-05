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
#     display_name: .venv2 (3.13.2.final.0)
#     language: python
#     name: python3
# ---

# %% [markdown]
# ## Installation and Set Up

# %%
# %pip install chemaxon rdkit pandas matplotlib scikit-learn seaborn boruta uvicorn joblib

# %%
# %env CHEMAXON_LICENSE_SERVER=https://license.chemaxon.com
# %env CHEMAXON_LICENSE_SERVER_KEY= <insert your license key here>

# %%
import chemaxon as cxn
import sklearn
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from chemaxon import *
# from chemaxon import io
# from chemaxon import calculations
# from chemaxon import pandasutil
# from chemaxon import fingerprints
# from chemaxon import standardizer
# from chemaxon import structurechecker
# from chemaxon import search
#from chemaxon.calculations import chemterms
from boruta import BorutaPy
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.metrics import mean_absolute_error, r2_score, mean_squared_error
from sklearn.model_selection import train_test_split
from scipy.stats import pearsonr
import joblib

print('Version:       ', cxn.__version__)
#print('CCL version:   ', cxn.ccl_version())
#print('CCL build date:', cxn.ccl_build_date())
print('License:       ', [lic for lic in cxn.licenses() if lic['product'] == 'Python API'])
print("Scikit-learn version:", sklearn.__version__)
print("Numpy version:", np.__version__)
print("Pandas version:", pd.__version__)
print("Seaborn version:", sns.__version__)


# %%
# Small helpers for loading local data files.
def load_text(rel_path):
    """Return a data file's contents as a string (e.g. XML configs, .mol files)."""
    with open(rel_path) as f:
        return f.read()


def load_pickle(rel_path):
    """joblib.load a pickle file."""
    return joblib.load(rel_path)


# %% [markdown]
# # Environment & resource helpers

# %% [markdown]
# Although you may have already provided your license through a different mechanism, Jupyter notebooks provide a convenient way to provide license keys
# Your license key is found on your Chemaxon account. I have already executed a cell setting my license key value to an environment variable, and won't repeat it for security reasons.
#
# Additionally in the above section, we have already performed the relevant pip installs and imports of various libraries.

# %%
# %env CHEMAXON_LICENSE_SERVER_KEY= <enter_your_license_key_here>

# %% [markdown]
# We define small loaders (`load_text`, `load_pickle`) once here and reuse them everywhere for reading the local data files.

# %% [markdown]
# ## Import and Clean Molecules

# %% [markdown]
# ### CSV Import, conversion to Chemaxon Molecule Type and Fragment Cleanup

# %%
# Load the preprocessed training CSV (SMILES already included)
df = pd.read_csv('compound_and_data_files/CHEMBL5763_training.csv')
# Display the first few rows
df[['CMP_CHEMBL_ID', 'SMILES', 'BIOACT_PCHEMBL_VALUE']]

# %%
#Convert SMILES to chemaxon mol object
df['cxn_mol_raw'] = df['SMILES'].apply(cxn.io.import_mol)
df.head()

# %% [markdown]
# A function to import directly from SDF is also available: open_for_import\
# You can find the details here: [https://apidocs.chemaxon.com/python_api/apidocs/chemaxon/io/importer.html](https://apidocs.chemaxon.com/python_api/apidocs/chemaxon/io/importer.html)\
# We've stuck with the above option for maximum flexibility

# %% [markdown]
# ### Helper functions — RDKit <-> ChemAxon interconversion

# %%
# RDKit <-> ChemAxon interconversion helpers (defined up front so they are available throughout).
from rdkit import Chem  # RDKit chemistry
import chemaxon
from rdkit.DataStructs.cDataStructs import SparseBitVect


def convert_rdkit_to_cxn(rdkit_molecule, intermediate_format="smiles", kekulization=True):
    """
    Converts an object in the rdkit molecule type to a chemaxon molecule type.
    Uses the intermediate format type specified.
    Assumes that the user desires the kekule structure, and not the aromatic representation, by default.
    Differences in output depending on the type of intermediate format used may exist.
    These differences can arise from limitations of certain intermediate filetypes to contain certain information.
    A less obvious difference may lie in the differences in default RDKit implementation of export functions - e.g. the SMILES seem to be aromatized, while other formats do not
    """
    match intermediate_format:
        case "smiles":
            cxn_molecule = chemaxon.io.import_mol(Chem.MolToSmiles(rdkit_molecule, kekuleSmiles=kekulization))
        case "mrv":
            cxn_molecule = chemaxon.io.import_mol(Chem.MolToMrvBlock(rdkit_molecule, kekulize=kekulization))
        case "molV2000":
            cxn_molecule = chemaxon.io.import_mol(Chem.MolToMolBlock(rdkit_molecule, kekulize=kekulization))
        case "molV3000":
            cxn_molecule = chemaxon.io.import_mol(Chem.MolToMolBlock(rdkit_molecule, forceV3000=True, kekulize=kekulization))
        case _:
            raise ValueError("Invalid intermediate_format: format must be smiles (default), mrv, molV2000 or molV3000.")

    return cxn_molecule


def convert_cxn_to_rdkit(cxn_molecule, intermediate_format="smiles", sanitization=True):
    """
    Converts an object in the cxn molecule type to an rdkit molecule type.
    Differences in output depending on the type of intermediate format used may exist.
    Defaults to sanitization=True as per the rdkit default options, but this can be modified.
    These differences can arise from limitations of certain intermediate filetypes to contain certain information.
    A less obvious difference may lie in the differences in default cxn implementation of export functions - e.g. the SMILES seem to be aromatized, while other formats do not
    """
    match intermediate_format:
        case "smiles":
            rdkit_molecule = Chem.MolFromSmiles(chemaxon.io.export_mol(cxn_molecule, 'smiles:u'), sanitize=sanitization)
        case "mrv":
            rdkit_molecule = Chem.MolFromMrvBlock(chemaxon.io.export_mol(cxn_molecule, 'mrv'), sanitize=sanitization)
        case "molV2000":
            rdkit_molecule = Chem.MolFromMolBlock(chemaxon.io.export_mol(cxn_molecule, 'mol:V2'), sanitize=sanitization)
        case "molV3000":
            rdkit_molecule = Chem.MolFromMolBlock(chemaxon.io.export_mol(cxn_molecule, 'mol:V3'), sanitize=sanitization)
        case _:
            raise ValueError("Invalid intermediate_format: format must be smiles (default), mrv, molV2000 or molV3000.")

    return rdkit_molecule


def to_sparse_vector(fp) -> SparseBitVect:
    """
    Converts from a chemaxon fingerprint to a SparseBitVect(or).
    These are used by certain RDKit fingerprint comparison functions, such as 'BulkTanimotoSimilarity'.
    """
    bv = SparseBitVect(len(fp.array) * 64)
    for int_index, int_val in enumerate(fp.array):
        for bit_index in range(64):
            if int_val & (1 << bit_index):
                bv.SetBit(int_index * 64 + bit_index)
    return bv


# %% [markdown]
# ### What do we do if there are multiple fragments present?
# It's a good idea to investigate the structures further in this case.\
# A common move is to remove all structures but the largest, which can be achieved with the "Remove Fragment" Standardizer action.\
# This isn't always perfect though. What if you have large salts? What if the entry genuinely represents a mixture?\
# If you have a fixed library of salts and solvates, you can use them in the "Strip Salts" or "Strip Solvates" actions.\
# This allows you to remove predetermined structures, and anything that still has multiple fragments afterwards may warrant further investigation

# %%
#Calculate total number of fragments. Let's see if we have any salts/solvates/other additional fragments present
df['fragment_count'] = df['cxn_mol_raw'].apply(lambda x: cxn.calculations.evaluate(x, "fragmentCount()"))
df.head()

# %%
#let's look at the multi_fragment structures
multi_fragment_df = df[df['fragment_count'].astype(int)!=1]
if multi_fragment_df.empty:
    print('Looks like everything is a single structure')
else:
    multi_fragment_df.head()

# %% [markdown]
# ### Standardizer - single structure walkthrough

# %% [markdown]
# Simple standardizer configurations can be specified with an input string.\
# For those with knowledge of the tool, this [configuration strings](https://docs.chemaxon.com/latest/standardizer_actions.html) can easily be input.

# %%
# Import a molecule. Notice the explicit Histidine represented as a contracted S-Group
structure_clean_with_strings_input = cxn.io.import_mol(load_text('compound_and_data_files/examorelin_s_group.mol'))
structure_clean_with_strings_input

# %%
# Let's get that structure cleaned up
structure_clean_with_strings_output = cxn.standardizer.Standardizer('ungroupsgroups').standardize(structure_clean_with_strings_input)
structure_clean_with_strings_output

# %% [markdown]
# For those who are newer to the tool, it's often easier to use the simple UI provided by the desktop tool to configure the actions.\
# This is also true as action sets become more complex, as the interface can make it easier to think through.\
# Otherwise, standardizer configuration XMLs can be visually created visually with the desktop UI.
#
# ![Screenshot of the Standardizer configuration pane](configurations/standardizer_config_screen.png)

# %%
# Configure a standardizer using a previously saved XML
standardizer_recommended_configuration = load_text('configurations/standardizers.xml')
standardizer_recommended = cxn.standardizer.Standardizer(standardizer_recommended_configuration)

# %%
# Import a molecule. Notice the explicit hydrogen, uncharged nitro group and methyl S-group
example_molecule_preclean = cxn.io.import_mol(load_text('compound_and_data_files/example_pre_clean.mol'))
example_molecule_preclean

# %%
# Let's clean it up
example_molecule_post_clean = standardizer_recommended.standardize(example_molecule_preclean)
example_molecule_post_clean

# %% [markdown]
# ### Standardizer - Standardize the whole data set

# %%
# Standardize the whole set
df['cxn_mol_standardized'] = df['cxn_mol_raw'].apply(lambda x: standardizer_recommended.standardize(x))
df.head()

# %% [markdown]
# #### Structure checking works much the same way
# There's a desktop UI for config, and you can use action strings or XMLs to define the configuration.\
# Just like with Standardizer, actions are configured in a predetermined order.\
# A notable difference is that for Structure Checkers, you can configure them to try to fix structures (in many cases), or simply identify the issue and take no action (in all cases).

# %%
bad_wedge = cxn.io.import_mol(load_text('compound_and_data_files/bad_wedge.mol'))
bad_wedge

# %%
wedge_flipper = cxn.structurechecker.StructureChecker('nonstereowedgebond->flipwedge')
fixed_wedge = wedge_flipper.fix(bad_wedge).fixed_mol
fixed_wedge

# %% [markdown]
# Sometimes, we want to inspect the process and results of a structure checker, and not simply return the (hopefully) fixed structure.\
# This is especially true when we are running a series of checks, or are trying to get a general feeling for the cleanliness of our input data set.\
# As well as the "fixed_mol" attribute, we can check for some other information

# %%
wedge_flipper_inspection = wedge_flipper.fix(bad_wedge)
wedge_flipper_inspection.is_fix_successful

# %% [markdown]
# Running the fix returns StructureCheckerBatchResult object\
# This objects attribute "results" is a list of StructureCheckerResult objects\
# You can check individual checker descriptions, atom indices, etc. here

# %%
#For more detailed results, check the StructureChecker object, rather than the Result object
bond_issue_inspection = wedge_flipper.check(bad_wedge)
print(str(bond_issue_inspection.results[0].description) + "\nOffending bond index " + str(bond_issue_inspection.results[0].bond_indices))


# %% [markdown]
# Let's return to fixing all of our input structures

# %%
# Structure checking
# Create a function to get both pieces of information that we need out
def structure_checker_structure_and_status_columns(
        raw_structure: cxn.Molecule,
        config_rel_path: str = "configurations/structure_checkers.xml"):
    str_chk_config = load_text(config_rel_path)
    fixed = cxn.structurechecker.StructureChecker(str_chk_config).fix(raw_structure)
    return (fixed.fixed_mol, fixed.is_fix_successful)


# %%
df[['cxn_mol_standardized_and_checked', 'check_fix_result']] = df['cxn_mol_standardized'].apply(
    lambda x: pd.Series(structure_checker_structure_and_status_columns(x))
)
df.head()

# %%
structures_failed_to_fix = df[~df['check_fix_result']]
if structures_failed_to_fix.empty:
    print('Looks like all of our structures were fixed succesfully')
else:
    structures_failed_to_fix.head()


# %%
#Let's overwrite the SMILES with a new SMILES output from our cleaned up structures, and remove some unnecessary columns
df['SMILES'] = df['cxn_mol_standardized_and_checked'].apply(lambda x: cxn.io.export_mol(x, 'SMILES'))
df.drop(columns = ['cxn_mol_raw', 'cxn_mol_standardized', 'fragment_count', 'check_fix_result'], inplace=True)
df.rename(columns={'cxn_mol_standardized_and_checked': 'cxn_mol'}, inplace=True)
df.head()

# %% [markdown]
# ## Microspecies determination
# The exact protonation state and tautomerization of a structure in a biological setting is important, and many descriptor generations are sensitive to them.\
# So, let's add one additional column featuring the major microspecies of the major tautomeric form.

# %%
df['cxn_mol_majms'] = df['cxn_mol'].apply(lambda x: cxn.calculations.major_microspecies(mol=x, ph = 7.4, take_major_tautomeric_form=True, keep_explicit_hydrogens=False))
df.head()


# %% [markdown]
# We keep the original form too as we find that not all predictors accept charged forms.
#
# Of course, it may be worthwhile investigating the full range of tautomers and microspecies available, especially if you find outliers, using the dedicated predictors separately.

# %% [markdown]
# ## Clean Data

# %%
#A quick visualization is always a good idea
sns.histplot(df, x=df['BIOACT_PCHEMBL_VALUE'].astype(float), kde=True, stat="count")
plt.show()

# %% [markdown]
# We should also check if there are duplicate structures present in the dataset.\
# Using a INCHIKey based text comparison is a valid approach, however we'll use the Chemaxon duplicate search since it enables greater control over tautomerization, target standardization and canonicalization.

# %%
self_duplicate_example = df['cxn_mol_majms'][42]

# Here, the default standardization of general aromatization will be applied to all compounds being compared
dup_search = cxn.search.MoleculeSearch(search_type=cxn.search.SearchType.DUPLICATE)

hit_list = dup_search.find_in_list(self_duplicate_example, df['cxn_mol_majms'].to_list())
indices = [i for i, val in enumerate(hit_list) if val]
df.iloc[indices]

# %% [markdown]
# If we run the above on a combinatorial search, it can take a few minutes.\
# In most situations, e.g. a database of compounds, you can run the searches in smaller chunks, i.e. as the compounds are being added.
#
# An alternative approach is a quick fingerprint screen, as shown below.\
# The Chemaxon Chemically Hashed Fingerprint, or CFP for short, is a proprietary fingerprint that is optimized for substructure and duplicate searches.

# %%
#A fingerprint comparison is a quicker process, however it sometimes turns up false positive matches
fp_search_copy = df.copy()
fp_search_copy['search_fps'] = fp_search_copy['cxn_mol_majms'].apply(lambda x: cxn.fingerprints.cfp(cxn.standardizer.Standardizer('aromatize').standardize(x)).bits)
potential_duplicate_mask = fp_search_copy.duplicated(subset=['search_fps'])
print("There may be as many as " + str(len(fp_search_copy[potential_duplicate_mask])) + " duplicates in this data set")

# %% [markdown]
# Since our search above did turn up some positives, we need to run a targeted, intensive search.\
# The full run of duplicate searches is also included for your reference below.\
# If you choose to run it, note that it may take a few minutes.\
# However, it is more exhaustive and less likely to run into false positive hits than the method we used above

# %%
# Run a proper duplicate search, but only on the previously identified positives.
dup_search_df_copy = df.copy()
copy_of_structures = dup_search_df_copy['cxn_mol_majms']
#set to -1 as we expect a self duplicate
dup_search_df_copy['num_duplicates'] = -1
for structure in copy_of_structures[potential_duplicate_mask]:
    hit_list = dup_search.find_in_list(structure, dup_search_df_copy['cxn_mol_majms'].to_list())
    indices = [i for i,val in enumerate(hit_list) if val]
    dup_search_df_copy.loc[indices, 'num_duplicates'] += 1

if max(dup_search_df_copy['num_duplicates']) == 0:
    print("We have no duplicates in the set")
else:
    print("We have duplicates in the set!")


# %% [markdown]
# ### Optional: Full duplicate search of all structures

# %%
# A full run of duplicates searches using the duplicate search type
dup_search_df_copy = df.copy()
copy_of_structures = dup_search_df_copy['cxn_mol_majms']
#set to -1 as we expect a self duplicate
dup_search_df_copy['num_duplicates'] = -1
for structure in copy_of_structures:
    hit_list = dup_search.find_in_list(structure, dup_search_df_copy['cxn_mol_majms'].to_list())
    indices = [i for i,val in enumerate(hit_list) if val]
    dup_search_df_copy.loc[indices, 'num_duplicates'] += 1

print(str(max(dup_search_df_copy['num_duplicates'])) + ", " + str(min(dup_search_df_copy['num_duplicates'])))

# %% [markdown]
# ## Out of Scope: Further Analysis of the data and compounds

# %% [markdown]
# Additional investigation is of course possible. This includes, but is not limited to:
# * Closer looks at the data distributions and their dynamic range
#     * Normalizing data if needed
# * Clustering the structures or reviewing chemical space
# * Scaffold analysis, from the clustering or otherwise
# * Reviewing duplicates in the data set, and adjusting, averaging or deleting multiply present values\
# \
# Given the time available to use today, this is unfortunately out of scope

# %% [markdown]
# ## Generate Descriptors

# %% [markdown]
# In the code below, we are defining a large number of descriptors.\
# There is an intentional overlap between Chemaxon, RDKit and other descriptors, including intentionally redundant ones\
# In practice, it is probably wiser to not intentionally calculate e.g. the number of aromatic rings multiple times\
# However, this is included here to demonstrate the interoperability of the two\
# \
# We'll take a look at reducing descriptor counts and removing redundancy later on

# %%
from descriptors import rdkit_descriptor_names

# %% [markdown]
# ### Dedicated calculation functions

# %% [markdown]
# A large number of our most popular predictors have their own dedicated functions that provide a structured output.\
# The below pka calculations has 9 parameters plus the structure available; we stick to the defaults here

# %%
viagra_pka = cxn.calculations.pka(cxn.io.import_mol('viagra'))
#Process the result for a cleaner output.
viagra_pka_df = pd.DataFrame(columns = ['atom_index', 'pka_type', 'pka_value'])
acidic_or_basic = ["Acidic", "Basic"]
for x in viagra_pka.pka_values:
    next_entry = {}
    next_entry['atom_index'] = x.atom_index
    next_entry['pka_type']=acidic_or_basic[x.pka_type]
    next_entry['pka_value']=x.value
    next_entry_df = pd.DataFrame([next_entry])
    viagra_pka_df = pd.concat([viagra_pka_df, next_entry_df])
viagra_pka_df.head()

# %% [markdown]
# Our ["Chemical Terms"](https://docs.chemaxon.com/latest/chemical-terms_functions-by-categories.html) provides a convenient catch-all for all other calculations.\
# Indeed, in this work we found it was more convenient to use Chemical Terms for bulk generation of many values, especially when a simple numeric result rather than a complex output was desired.

# %%
from rdkit import Chem
from rdkit.ML.Descriptors import MoleculeDescriptors
#from descriptors import rdkit_descriptor_names

rdkit_calculator = MoleculeDescriptors.MolecularDescriptorCalculator(rdkit_descriptor_names)

def compute_properties_from_mol(cxn_molecule):
    
    mol = cxn_molecule
    rdkit_mol = convert_cxn_to_rdkit(cxn_molecule=cxn_molecule)

    #CHEMTERM helper function
    def c(expr):
        return cxn.calculations.evaluate(mol, expr)
    
    # CHEMAXON DESCRIPTORS ---------------

    # Generate logD at multiple pH values
    logd_ph_values = [3.4, 5.4, 7.4, 9.4, 11.4]
    logd_dict = {f'LOGD[{pH}]': cxn.calculations.logd(mol, ph=pH) for pH in logd_ph_values}

    # Generate atomcounts for chosen atoms
    atomcount_values = [1, 6, 7, 8, 9, 15, 16, 17]
    atomcount_dict = {f'ATOMCOUNT({No})': c(f"atomCount('{No}')") for No in atomcount_values}

    # Generate ECFP fingerprint bit array
    fp = cxn.fingerprints.ecfp(mol, 4, 1024)
    binary_string = ''.join(fp.to_binary_string())
    bit_array = list(map(int, binary_string))

    # RDKIT DESCRIPTORS ---------------

    # RDkit descriptor calculation
    rdkit_descriptor_values = rdkit_calculator.CalcDescriptors(rdkit_mol)
    rdkit_dict = dict(zip(rdkit_descriptor_names, rdkit_descriptor_values))

    # DESCRIPTOR GENERATION ---------------

    return pd.Series({
        **logd_dict,
        **atomcount_dict,
        **{f'FP_{i}': bit for i, bit in enumerate(bit_array)},
        'LOGP': cxn.calculations.logp(mol),
        'PSA': c("psa()"),
        'apKamin': c("min(apKa('1'),14)"),
        'apKamax': c("max(bpKa('1'),0)"),
        'hlb': c("hlb()"),
        'mass': c("mass()"),
        'maxcharge': c("max(charge())"),
        'mincharge': c("min(charge())"),
        'donorCount': c("donorCount()"),
        'donorSiteCount': c("donorSiteCount()"),
        'acceptorCount': c("acceptorCount()"),
        'acceptorSiteCount': c("acceptorSiteCount()"),
        'hmoPiEnergy()': c("hmoPiEnergy()"),
        'mpo': c("mpo()"),
        'logS': c("logS()"),
        'aliphaticAtomCount': c("aliphaticAtomCount"),
        'aliphaticRingCount': c("aliphaticRingCount"),
        'aromaticAtomCount': c("aromaticAtomCount"),
        'aromaticRingCount': c("aromaticRingCount"),
        'atomCount': c("atomCount"),
        'balabanIndex': c("balabanIndex"),
        'bondCount': c("bondCount"),
        'carboRingCount': c("carboRingCount"),
        'chainAtomCount': c("chainAtomCount"),
        'cyclomaticNumber': c("cyclomaticNumber"),
        'fsp3': c("fsp3"),
        'fusedAliphaticRingCount': c("fusedAliphaticRingCount"),
        'fusedAromaticRingCount': c("fusedAromaticRingCount"),
        'fusedRingCount': c("fusedRingCount"),
        'hararyIndex': c("hararyIndex"),
        'heteroaliphaticRingCount': c("heteroaliphaticRingCount"),
        'heteroaromaticRingCount': c("heteroaromaticRingCount"),
        'heteroRingCount': c("heteroRingCount"),
        'hyperWienerIndex': c("hyperWienerIndex"),
        'largestConjugatedSystemSize': c("largestConjugatedSystemSize"),
        'largestRingSize': c("largestRingSize"),
        'largestRingSystemSize': c("largestRingSystemSize"),
        'plattIndex': c("plattIndex"),
        'randicIndex': c("randicIndex"),
        'ringAtomCount': c("ringAtomCount"),
        'ringBondCount': c("ringBondCount"),
        'ringSystemCount': c("ringSystemCount"),
        'rotatableBondCount': c("rotatableBondCount"),
        'smallestRingSize': c("smallestRingSize"),
        'smallestRingSystemSize': c("smallestRingSystemSize"),
        'szegedIndex': c("szegedIndex"),
        'wienerIndex': c("wienerIndex"),
        'wienerPolarity': c("wienerPolarity"),
        **rdkit_dict
    })


# %% [markdown]
# Below, we are calculating descriptors based on the SMILES column. We've chosen to do this for simplicity, and because this is commonplace in public examples.\
# The Chemaxon calculators can of course generate values based on the Chemaxon molecule object. \
# If that's of more interest, we have provided some helper functions at the top of the notebook that facilitate quick interconversion between Chemaxon and RDKit types using an intermediate cheminformatic format (the type of which you can choose in the function parameters)

# %%
#Calculate descriptors
df_props = df['cxn_mol_majms'].apply(compute_properties_from_mol)

# Combine results with original DataFrame
df = df.drop(columns=[col for col in df_props.columns if col in df.columns], errors='ignore')
df = pd.concat([df.drop(columns=[col for col in df_props.columns if col in df.columns], errors='ignore'),
                df_props], axis=1)

pd.set_option('display.max_columns', None)

df

# %% [markdown]
# ### Optional: read (or write) pre-calculated descriptors from a CSV to save time

# %% [markdown]
# The above cell can take a while to run, so we've saved the output into a CSV.\
# To save time, you might wish to load it in instead of re-running the cell!\
# This will replace the cxn_mol column with a Chemaxon Extended Smiles string

# %%
# Saving DF as CSV for saving time
df.to_csv('compound_and_data_files/calc_descriptors.csv', index=False)


# %%
df = pd.read_csv('compound_and_data_files/calc_descriptors.csv')
df['cxn_mol'] = df['SMILES'].apply(cxn.io.import_mol)
df

# %% [markdown]
# At this point, it might be a good idea to take a look at the chemical space that our data set occupies, especially if you have a dedicated test set to compare it to.\
# At the very least, comparing future compounds you wish to predict against the training set gives an idea of whether they fall within the applicability domain.\
# A quick, convenient proxy is to take a new compound and retrieve its highest similarity score against the training set, as shown below.\
# We show two views: an **ECFP / Tanimoto** view (substructure-style similarity) and a **pharmacophore-fingerprint** view (binding-feature similarity).

# %%
aspirin = cxn.io.import_mol('aspirin')

# 1) ECFP (circular) similarity
ecfp_query = cxn.fingerprints.ecfp(aspirin, 4, 1024)
similarities = df[['cxn_mol']].copy()
similarities['ecfp_similarity'] = similarities['cxn_mol'].apply(
    lambda m: cxn.fingerprints.tanimoto(ecfp_query, cxn.fingerprints.ecfp(m, 4, 1024))
)

# 2) Pharmacophore-fingerprint similarity (binding-feature based)
pfp_query = cxn.fingerprints.pharmacophore_fp(aspirin)
similarities['pharmacophore_similarity'] = similarities['cxn_mol'].apply(
    lambda m: cxn.fingerprints.float_vector_tanimoto(pfp_query, cxn.fingerprints.pharmacophore_fp(m))
)

closest_matches = similarities.nlargest(3, 'ecfp_similarity')
closest_matches.head()


# %% [markdown]
# ## Reduce Descriptor Count

# %% [markdown]
# We need to remove some of the highly correlated features

# %%
def remove_highly_correlated_features(df, threshold=0.95):
    corr_matrix = df.corr().abs()
    upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))
    to_drop = [column for column in upper.columns if any(upper[column] > threshold)]
    return df.drop(columns=to_drop)

drop_cols = ['CMP_CHEMBL_ID', 'ID', 'Molecule', 'SMILES', 'BIOACT_PCHEMBL_VALUE', 'cxn_mol', 'cxn_mol_majms']
x = df.drop(columns=[col for col in drop_cols if col in df.columns])
y = df['BIOACT_PCHEMBL_VALUE']

x = remove_highly_correlated_features(x)
x

# %% [markdown]
# It's also a good idea to remove descriptors that don't make a difference to the model. The Boruta algorithm is a good way to do this without having to set a variance cut off.\
# We'll do this after splitting into test and train to prevent any data leakage.

# %%
# Train-test split
x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.2, random_state=832)

# %%
# Use a simple RF Regressor for speedy descriptor selection
boruta_regressor = RandomForestRegressor(n_jobs=1, max_depth=5)

# Boruta keeps only features that consistently beat randomised "shadow" copies of themselves
feature_selector = BorutaPy(boruta_regressor, n_estimators='auto', max_iter=20, random_state=67, perc=75)
feature_selector.fit(np.array(x_train), np.array(y_train))

# %%
#This produces a filtered set of features
feature_selector.support_

# %%
x_train_filtered = feature_selector.transform(x_train, return_df=True)
x_train_filtered

# %% [markdown]
# ### Optional: Write/Read to/from pickle object to save time

# %% [markdown]
# All of the above can take some time, so we've again made it easy to save and load.

# %%
# Write to pickle
joblib.dump(x_train_filtered.columns.tolist(), "filtered_features.pkl")

# %%
# Read from pickle
list_of_features = load_pickle('filtered_features.pkl')
x_train_filtered = x_train[list_of_features]
x_train_filtered

# %% [markdown]
# # Filter test set

# %%
x_test_filtered = x_test.filter(x_train_filtered.columns.tolist())
x_test_filtered

# %% [markdown]
# ## Train Model(s)

# %%
### 1. RANDOM FOREST
rf = RandomForestRegressor(
    n_estimators=100,
    max_features='sqrt',
    max_depth=25,
    min_samples_leaf=2,
    bootstrap=True,
    random_state=475
)
rf.fit(x_train_filtered, y_train)
y_pred_rf = rf.predict(x_test_filtered)

mae_rf = mean_absolute_error(y_test, y_pred_rf)
rmse_rf = np.sqrt(mean_squared_error(y_test, y_pred_rf))
r2_rf = r2_score(y_test, y_pred_rf)
pcc_rf, _ = pearsonr(y_test, y_pred_rf)

### 2. GRADIENT BOOSTING
gb = GradientBoostingRegressor(
    n_estimators=300,
    max_depth=25,
    min_samples_leaf=5,
    max_leaf_nodes=20,
    subsample=0.7,
    learning_rate=0.05,
    loss='squared_error',
    random_state=475
)

gb.fit(x_train_filtered, y_train)
y_pred_gb = gb.predict(x_test_filtered)

mae_gb = mean_absolute_error(y_test, y_pred_gb)
rmse_gb = np.sqrt(mean_squared_error(y_test, y_pred_gb))
r2_gb = r2_score(y_test, y_pred_gb)
pcc_gb, _ = pearsonr(y_test, y_pred_gb)

### Print comparison
print("Random Forest:")
print(f"  MAE : {mae_rf:.3f}")
print(f"  RMSE: {rmse_rf:.3f}")
print(f"  R²  : {r2_rf:.3f}")
print(f"  PCC : {pcc_rf:.3f}")

print("\nGradient Boosting:")
print(f"  MAE : {mae_gb:.3f}")
print(f"  RMSE: {rmse_gb:.3f}")
print(f"  R²  : {r2_gb:.3f}")
print(f"  PCC : {pcc_gb:.3f}")

# %% [markdown]
# ## Evaluate Model(s)

# %% [markdown]
# Predicted vs. experimental scatter plot

# %%
# RF model
sns.scatterplot(x=y_test, y=y_pred_rf)
lims = [
    min(y_test.min(), y_pred_rf.min()),
    max(y_test.max(), y_pred_rf.max()),
]
plt.plot(lims, lims, color="red", linestyle="--", linewidth=1)
plt.xlim(lims)
plt.ylim(lims)
plt.xlabel('Experimental')
plt.ylabel('Predicted')
plt.title('RF model')
plt.show()

# GB model
sns.scatterplot(x=y_test, y=y_pred_gb)
lims = [
    min(y_test.min(), y_pred_gb.min()),
    max(y_test.max(), y_pred_gb.max()),
]
plt.plot(lims, lims, color="red", linestyle="--", linewidth=1)
plt.xlim(lims)
plt.ylim(lims)
plt.xlabel('Experimental')
plt.ylabel('Predicted')
plt.title('GB model')
plt.show()

# %% [markdown]
# In the next cell we select whichever model scored the higher **R²** on the test set (done programmatically). In this run that is typically the Gradient Boosting model.\
# Note that comparing two models on a single train/test split like this is **not** a statistically rigorous selection procedure — n-fold cross-validation with a proper statistical comparison is the correct approach (see *Out of scope: Refine Model(s)* below).

# %%
# Choose model based on R² score
if r2_rf >= r2_gb:
    best_model = rf
    print("\nSelected model: Random Forest")
else:
    best_model = gb
    print("\nSelected model: Gradient Boosting")

# Save the model
joblib.dump(best_model, 'best_model.pkl')

# %% [markdown]
# ## Out of scope: Refine Model(s)

# %% [markdown]
# There are a number of additional steps that can be taken here, which include
# * Testing further model types and architecture
# * Hyperparameter tuning
# * More thorough approaches like n-fold cross validation and statistical comparison of different models generated

# %% [markdown]
# ## Predict Property(ies) with Model(s)

# %%
# Import the held-out prediction set
df_predict = pd.read_csv('compound_and_data_files/CHEMBL5763_pred.csv')
df_predict.head()

# %%
# Compute properties per SMILES
df_predict['cxn_mol'] = df_predict['SMILES'].apply(lambda x: cxn.io.import_mol(x))
df_predict['cxn_mol_majms'] = df_predict['cxn_mol'].apply(lambda x: cxn.calculations.major_microspecies(mol=x, ph = 7.4, take_major_tautomeric_form=True, keep_explicit_hydrogens=False))
df_props_predict = df_predict['cxn_mol_majms'].apply(compute_properties_from_mol)

# Combine results with original DataFrame
df_predict = df_predict.drop(columns=[col for col in df_props_predict.columns if col in df_predict.columns], errors='ignore')
df_predict = pd.concat([df_predict, df_props_predict], axis=1)

pd.set_option('display.max_columns', None)
df_predict.head()

# %%
# Load the trained model from file
loaded_model = load_pickle('best_model.pkl')

# %%
#Filter the features 
df_features = df_predict[list_of_features]
df_features.head()

# %%
# Predict
y_pred = loaded_model.predict(df_features)

# Add predictions to the DataFrame
df_result = df_predict[['SMILES', 'BIOACT_PCHEMBL_VALUE']].copy()
df_result['Predicted_BIOACT_PCHEMBL_VALUE'] = y_pred

df_result.head()

# %%
y_exp = df_result['BIOACT_PCHEMBL_VALUE']
y_pred = df_result['Predicted_BIOACT_PCHEMBL_VALUE']

mae = mean_absolute_error(y_exp, y_pred)
rmse = np.sqrt(mean_squared_error(y_exp, y_pred))
r2 = r2_score(y_exp, y_pred)
pcc, _ = pearsonr(y_exp, y_pred)

print(f"\nModel performance on the new set of compounds")
print(f"  MAE : {mae:.3f}")
print(f"  RMSE: {rmse:.3f}")
print(f"  R²  : {r2:.3f}")
print(f"  PCC : {pcc:.3f}")

# %%
sns.scatterplot(x=y_exp, y=y_pred)
lims = [
    min(y_exp.min(), y_pred.min()),
    max(y_exp.max(), y_pred.max()),
]
plt.plot(lims, lims, color="red", linestyle="--", linewidth=1)
plt.xlim(lims)
plt.ylim(lims)
plt.xlabel('Experimental')
plt.ylabel('Predicted')
plt.title('Predicted vs experimental BIOACT_PCHEMBL_VALUE')
plt.show()

# %% [markdown]
# ## How can we distribute this to users?

# %% [markdown]
# Plugins can easily be created by exposing a rest API.\
# If we take a look into the ML_model_plugin.py file, we can see an example a plugin that is compatible with our Design Hub design management platform.

# %% [markdown]
# ## Additional Feature Spotlight
#
# Key features of the Python API that didn't fit into our session today

# %% [markdown]
# ### Built-in ADMET predictors
#
# Ready-made models you would otherwise build from scratch — useful as features *or* as baselines for the model we just trained: hERG (classification + activity), blood-brain-barrier penetration, and CNS MPO.

# %%
demo_mol = cxn.io.import_mol('warfarin')
print('hERG class   :', cxn.calculations.herg_classification(demo_mol))
print('hERG activity:', cxn.calculations.herg_activity(demo_mol))
print('BBB          :', cxn.calculations.bbb(demo_mol))
print('CNS MPO      :', cxn.calculations.cns_mpo(demo_mol))

# %% [markdown]
# ### Combinatorial library enumeration (Reactor)
#
# Beyond a static dataset, `chemaxon.reactor` enumerates virtual libraries from a reaction plus reactant lists — a natural "where to go next": generate new candidates, featurize them with the function above, and score them with the trained model. The shape is roughly:
#
# ```python
# # Amide coupling: carboxylic acid + amine -> amide
# reaction = cxn.import_mol('[C:1](=O)O.[N:2]>>[C:1](=O)[N:2]')
# acids  = [cxn.import_mol(s) for s in ['CC(=O)O', 'OC(=O)c1ccccc1']]
# amines = [cxn.import_mol(s) for s in ['CCN', 'NCc1ccccc1']]
#
# result = cxn.reactor.parallel(reaction, [acids, amines], mode=cxn.reactor.Mode.COMBINATORIAL)
# ```
#
# Reactor's `ReactorOptions` / `ReactionRules` also control reactivity, selectivity, tolerances and atom-mapping.
# Full details [here](https://github.com/ChemAxon/python-examples/blob/main/jupyter/09_reactor.ipynb)

# %% [markdown]
# ## Helper Functions

# %%
#Make use of the pandasutils to display mols as images in dfs
from IPython.display import display, HTML
# df.to_html('web_view.html', escape=False, formatters=dict(cxn_mol=cxn.mol_to_svg_formatter))
display(HTML(df.to_html('web_view.html', escape=False, formatters=dict(cxn_mol=cxn.pandasutil.mol_to_svg_formatter))))

# %% [markdown]
# As well as the documentation pages and API reference, we have simple Python files and Jupyter notebooks published on our GitHub to help you gain an understanding of the API

# %% [markdown]
# [Documentation page](https://docs.chemaxon.com/latest/python-api_index.html)\
# [API reference](https://apidocs.chemaxon.com/python_api/apidocs/chemaxon.html)\
# [Public examples page](https://github.com/ChemAxon/python-examples/tree/main)

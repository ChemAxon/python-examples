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
# ### **Reactor module**

# %%
#Helper function to show the reactor results

def display_result(result: list):
    for r in range(len(result)):
        print(r + 1, '. result', sep='')
        for i in range(len(result[r].product_sets)):
            print('  ', i + 1, '. product set:', sep='')
            for j in range(len(result[r].product_sets[i].products)):
                print('    ', j + 1, '. product', sep='')
                display(result[r].product_sets[i].products[j])


# %%
from chemaxon.io import import_mol, export_mol
from chemaxon.reactor import react_single_input, react_parallel, Mode, ReactorOptions, ResultType, ReactionRules

reactants1 = [import_mol('CC(O)=O'), import_mol('CCC(O)=O')]
reactants2 = [import_mol('CO'), import_mol('CCCCO')]
reaction = import_mol('[H:5][O:4][CH:1]=[O:2].[H:7][O:6][CH3:8]>>[CH3:8][O:6][CH:1]=[O:2]')
display(reaction)

# %%
result = react_single_input(reaction, [reactants1[0], reactants2[0]])
display_result([result])

# %%
result = react_parallel(Mode.SEQUENTIAL, reaction, [reactants1, reactants2])
display_result(result)

# %%
result = react_parallel(Mode.COMBINATORIAL, reaction, [reactants1, reactants2])
display_result(result)

# %%
all_reactants = [[reactants1[0], reactants1[1]], [reactants1[0], reactants2[1]]]
result = react_parallel(Mode.SEQUENTIAL, reaction, all_reactants)
print('Results without reaction rules')
display_result(result)

print('Results with reaction rules')
rules = ReactionRules()
rules.reactivity = 'match(reactant(1), alcohol) or match(reactant(1), phenol)'
result = react_parallel(Mode.SEQUENTIAL, reaction, all_reactants, rules=rules)
display_result(result)

# %%
options = ReactorOptions()
options.reverse = True
options.result_type = ResultType.REACTION
reverse_reactants = [[import_mol('COC(C)=O'), import_mol('CCCCOC(C)=O')]]
result = react_parallel(Mode.SEQUENTIAL, reaction, reverse_reactants, options=options)
display_result(result)

# %%
reaction = import_mol('[H:5][O:4][CH:1]=[O:2].[H:7][O:6][CH3:8]>>[CH3:8][O:6][CH:1]=[O:2]')
all_reactant = [[import_mol('OC(=O)CCCCC(CCC(O)=O)CC(O)=O')], [import_mol('CCO')]]
options = ReactorOptions()
print('Ratio 1:1')
options.ratio = [1, 1]
result = react_parallel(Mode.SEQUENTIAL, reaction, all_reactant, options=options)
display_result(result)

print('Ratio 1:2')
options.ratio = [1, 2]
result = react_parallel(Mode.SEQUENTIAL, reaction, all_reactant, options=options)
display_result(result)

print('Ratio 1:3')
options.ratio = [1, 3]
result = react_parallel(Mode.SEQUENTIAL, reaction, all_reactant, options=options)
display_result(result)

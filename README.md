# Chemaxon Python API examples

This repository contains usage examples for [Chemaxon Python API](https://docs.chemaxon.com/display/docs/python-api_index.md).
For each notebook (ipynb) file the corresponding py file contain the pure python commands. This python file has a clear git history without generated cell outputs.

## Prerequisites

- Chemaxon license key (or license file). See the [documentation](https://docs.chemaxon.com/display/docs/python-api_installation.md#license-installation) for details.
- [JupyterLab](https://jupyter.org/install) for Jupyter notebook examples.

## Note

If you have any question, suggestion please feel free to contact us via
[Chemaxon Support Portal](https://chemaxon.freshdesk.com/a/tickets/new)

## Contributors 

- For contributors only: install [jupytext](https://jupytext.org/), which will enable synchronization of notebook and python (ipynb and py) files.
- If pairing newly created notebook and python files is needed: `ls jupyter/*.ipynb | xargs jupytext --set-formats ipynb,py:percent` will generate the python file and pair the notebook file with it.

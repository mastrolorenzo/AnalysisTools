## 1. Setup

If this is your first time setting up the workspace, running the `setup` script to check for an install any missing Python packages.

```bash
./setup
```

## 2. Usage

First, create a SkimTreeBuilder configuration file. The configuration file is simply a Python module which requires the following attributes to be provided:

* `input_pattern`
   The input pattern used to glob for sample files.  The input pattern must contain a "%s" which will be substituted by a sample's input token.
* `destination`
   The output directory containing the skimmed trees.
* `additional_branches`
    The list of branches to keep in addition to the default ones. If there are no additional branches you would like to save, just pass an empty list "[]".
* `input_tokens`
   A list of input tokens corresponding to their samples. These are substituted into the input pattern to glob for sample files.

Please see the file `config.py` for example settings and feel free to use it as a template.

When ready, run the SkimTreeBuilder using the `run` script which takes a configuration file path as its only argument. The script will then build skimmed trees for each sample whose input token is present in the configuration file and place them in the output directory specified by the same configuration file.

```bash
./run config.py
```

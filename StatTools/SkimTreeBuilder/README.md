## 1. Setup

If this is your first time setting up the workspace, running the `setup` script to check for an install any missing Python packages.

```bash
./setup
```

## 2. Usage

First, create a SkimTreeBuilder configuration file. The configuration file is simply a Python module which requires the following attributes to be defined:

* `input_pattern`
   The input pattern used to glob for sample files. The input pattern must contain a "%s" which will be substituted by a sample's input token.
* `destination`
   The output directory containing the skimmed trees.
* `keep_branches`
   The list of existing branches to keep in addition to the default subset of branches. Don't worry about accidentally including the same branch twice or a default branch in this list, as duplicates are handled correctly. If there are no additional branches to be kept, just pass an empty list "[]".
* `new_branches`
   The list of new branches to create for the skimmed trees. Each item in the list must be a function which takes as input a single event and returns the value for the new branch. Don't forget to add any branches that the function depends on to `keep_branches` if it isn't kept by default. The name of the function must follow the convention <branch\_name>\_<branch\_type>, where the first part of its name before the last underscore is used to set the name of the branch and the second part of its name after the last underscore is a single character which sets the type of the branch. The following basic types are supported:

   - `b` for 8 bit unsigned integer
   - `s` for 16 bit unsigned integer
   - `i` for 32 bit unsigned integer
   - `l` for 64 bit unsigned integer
   - `B` for 8 bit signed integer
   - `S` for 16 bit signed integer
   - `I` for 32 bit signed integer
   - `L` for 64 bit signed integer
   - `F` for 32 bit floating point
   - `D` for 64 bit floating point
   - `O` for boolean

   If the name of the new branch matches the name of an existing branch then the values of the existing branch will be overwritten and potentially cast to a different type. If there are no new branches to be created, just pass an empty list "[]".
* `samples`
   A dictionary keeping track of the samples for which skimmed trees are built. The keys are the names of the samples and are also used to name the files containing their skimmed trees.  of input tokens corresponding to their samples. The values are lists of arguments pertaining to each sample, which are in order:
   - The scale factors for the sample. This must be a function which takes as input a single event and modifies the value of its weight branch. Don't forget to add any branches that the function depends on to `keep_branches` if it isn't kept by default. If there are no scale factors to apply, set this to `None`.
   - The list of input tokens for the sample. These are substituted into the input pattern to glob for sample files. While there is usually just one input token for a sample, multiple input tokens can be provided if the sample is actually a combination of other samples.

For an example of basic usage, please see the module `config.py`.
For an example of scale factor, please see the module `scale_factors_HIG16044.py`.
For an example of new branches, please see the module `new_branches.py`.

Feel free to use these example configuration files as templates to be copied and modified as needed.

When ready, run the SkimTreeBuilder using the `run` script which takes a configuration file path as its only argument. The script will then build skimmed trees for each sample requested and place them in the specified output directory.

```bash
./run config.py
```

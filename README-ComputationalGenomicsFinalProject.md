## Final Project

### Overview

For our final project, we added features to a Python library for phylogenetic analysis called TreeTime. Our main improvements to TreeTime are the addition of rate variation and site-specific matrices for both joint and marginal ancestral sequence reconstruction. TreeTime supports a number of functions besides sequence reconstruction. For a fuller picture, here is how the main README.txt summarizes it:

> TreeTime provides routines for ancestral sequence reconstruction and inference of molecular-clock phylogenies, i.e., a tree where all branches are scaled such that the positions of terminal nodes correspond to their sampling times and internal nodes are placed at the most likely time of divergence. 
>
> To optimize the likelihood of time-scaled phylogenies, TreeTime uses an iterative approach that first infers ancestral sequences given the branch length of the tree, then optimizes the positions of unconstrained nodes on the time axis, and then repeats this cycle. The only topology optimization are (optional) resolution of polytomies in a way that is most (approximately) consistent with the sampling time constraints on the tree. The package is designed to be used as a stand-alone tool on the command-line or as a library used in larger phylogenetic analysis work-flows

When ancestral sequence reconstruction is used as a sub-step in complex analytical procedures, such as the expectation-maximization system mentioned above, rate variation and site-specific options can still be applied. However, we focused on improving standalone ancestral sequence reconstruction. 

### Installation and Prerequisites

How to install TreeTime, and the prerequisites for doing so, haven't changed. Here is the relevant section of the main README.txt:

> TreeTime is compatible with Python 2.7 upwards and is tested on 2.7, 3.5, and 3.6. It depends on several Python libraries:
> * numpy, scipy, pandas: for all kind of mathematical operations as matrix operations, numerical integration, interpolation, minimization, etc.
> * BioPython: for parsing multiple sequence alignments and all phylogenetic functionality
> * matplotlib: optional dependency for plotting
>
> You may install TreeTime and its dependencies by running
>```bash
>  pip install .
>```
> within this repository. You can also install TreeTime from PyPi via
>```bash
>  pip install phylo-treetime
>```
> You might need root privileges for system wide installation. Alternatively, you can simply use it TreeTime locally without installation. In this case, just download and unpack it, and then add the TreeTime folder to your $PYTHONPATH.

### Ancestral Sequence Reconstruction (Vanilla)

There are two steps to (standalone) ancestral sequence reconstruction in TreeTime. First, instantiate a TreeTime tree. Second, run ancestral sequence reconstruction on it. The 3 standard parameters for a tree are a topology, a multiple-sequence alignment (a MSA), and a general time-reversible model (a GTR). Topologies and MSAs are read in from files. Acceptable file formats for topologies are newick, nexus and phylip,  and for MSAs fasta and phylip. GTRs are specified with the name of a substitution model, e.g. "Jukes-Cantor" or "WAG"; the full set of options is found in the ```standard``` function in ```gtr.py```. This last choice is important because when a tree has ancestral sequence reconstruction or other phylogenetic techniques run on it, the substitution model that they use in their calculations is determined by the GTR attatched to the tree.

```python
from treetime import TreeAnc

# instatiate tree
myTree = TreeAnc(gtr="Jukes-Cantor", tree="data/mytree.nwk", aln="data/mytree.fasta")

# actually do reconstruction (marginal)
myTree.infer_ancestral_sequences(marginal=True)

# actually do reconstruction (joint)
myTree.infer_ancestral_sequences(marginal=False)
```

### Rate Variation
Rate distributions are always discrete in TreeTime. Trees store the rate distribution that applies to them as a list of ```(rate, log_probability_of_rate)``` tuples. Trees have an optional ```rates``` parameter to set the rate distribution that is used for them. ```rates``` expects an input in the standard form: a list of ```(rate, log_probability_of_rate)``` tuples. If not passed a rate-setting parameter, trees will default to the uniform distribution,```[(1, 0)]```. 

The ```ASRV``` class derives gamma distributions and presents them in a useful format for modelling rate variation. ```ASRV.calc_rates()``` returns an 8-rate approximation of the continuous gamma distribution parametrized by a given alpha value. The list returned by ```calc_rates()``` is correctly formatted for use in the ```rates``` parameter.

```python
import numpy as np
from treetime import TreeAnc
from treetime.asrv import ASRV

# manual rates
myTree = TreeAnc(rates=[(.5, np.log(1/3)), (1, np.log(1/3)), (1.5, np.log(1/3))], gtr="Jukes-Cantor", tree="data/mytree.nwk", aln="data/mytree.fasta")

# gamma-distributed rates from ASRV class
asrv = ASRV(alpha=.5)
rates = asrv.calc_rates()
myTree = TreeAnc(rates=rates, gtr="Jukes-Cantor", tree="data/mytree.nwk", aln="data/mytree.fasta")
```

If marginal or joint reconstruction is run on a tree that has a non-uniform rate distribution, we automatically switch to performing the reconstruction with an algorithm that takes multiple rates into account. The branch and bound algorithm is used under the hood for joint reconstruction with rate variation. 

### Site-Specific Transition Probabilities
A tree will have site-specific matrices if its GTR parameter is set one of ```"BE", "SECONDARY", "EX", "EHO"```. We think it fits to use the GTR parameter like this because site-specific matrices are just a very sophisticated substitution model in some sense. As would be expected from the way the GTR parameter is used in a regular tree, ancestral sequence reconstruction on a tree with site-specific matrices uses those matrixes to provide the substitution model (at their corresponding sites) that is used in the reconstruction calculations. ```"BE"``` and ```"EX"``` invoke the same matrices; ```"SECONDARY"``` and ```"EHO"``` are similarly paired. (Why there are "repeats" is explained below.) The two matrices of ```"BE"/"EX"``` encode solvent accessibility: whether a site is buried or exposed. The three matrices of ```"SECONDARY"/"EHO"``` encode protein secondary structure: whether a site is on the helix, extended away, or something else ("O" for "other"). Using these models naturally requires specifying the category to which each site in the input sequence belongs; this is done with an additional tree instantiation parameter, ```struct_propty```. ```struct_propty``` is an array (Python list) of the same length as the input sequence (i.e. as the MSA). The nth value in the array is a one-letter string corresponding to the catagory of nth site in the sequence: ```"B"``` for buried, ```"E"``` for exposed, ```"H"``` for helix, ```"E"``` again for extend, and ```"O"``` for other. The reason why there are "repeat" matrix classes is a trade-off. ```"BE"``` and ```"SECONDARY"``` are quite slow but the GTR property of a tree set to one of them retains its expected behavior. ```"EX"``` and ```"EHO"``` make opposite choice, breaking some parts of the inferface in exchange for speed. Site-specific matrices can be used in conjunction with rate variation. 

```python
from treetime import TreeAnc

# pretend input sequence is length 5
solvent_accessibility = ["B", "B", "E", "E," "B"] 

# create tree using buried-exposed matrices
myTree = TreeAnc(gtr="BE", struct_propty=solvent_accessibility, tree="data/mytree.nwk", aln="data/mytree.fasta")
```

  
### Changed Files
This is an exhaustive list of the files in TreeTime that we modified for our project.
* treetime/treeanc.py
* treetime/gtr_site_specific.py
* treetime/seq_utils.py
* treetime/asrv.py (created)
* treetime/run_test.py
* treetime/test_anc.py (created)
* test_treetime.py

### Extended Example
See the "TreeTime Demo" Jupyter notebook. For some reason, Jupyter notebooks don't seem to work unless TreeTime is fully installed; it isn't sufficient to just put the TreeTime folder on your PYTHONPATH.

The tests are good demos as well. In paticular, ```test_ancestral``` shows biological data being loaded from appropriately formatted files and then used in both joint and marginal reconstruction. 
 
### Related Tools
[PASTA](https://github.com/smirarab/pasta) for generating MSAs, [RAxML](https://github.com/stamatak/standard-RAxML) for generating tree toplogies and inferring alpha parameters, [DSSP](https://github.com/cmbi/hssp)/[Sable](http://sable.cchmc.org) for inferring protein structual information.
 

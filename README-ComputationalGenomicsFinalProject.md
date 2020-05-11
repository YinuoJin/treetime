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

There are two steps to (standalone) ancestral sequence reconstruction in TreeTime. First, instantiate a TreeTime tree. Second, run ancestral sequence reconstruction on it. The 3 standard parameters for a tree are a topology, a multiple-sequence alignment (a MSA), and a general time-reversible model (a GTR). Topologies and MSAs are read in from files. Acceptable file formats for topologies are newick, nexus and phylip,  and for MSAs fasta and phylip. GTRs are specified with the name of a substitution model, e.g. "Jukes-Cantor" or "WAG"; the full set of options is found in the ```standard``` function in ```gtr.py```. (A more advanced option is for users to create a custom GTR, which can be parametrized by an arbitrary substitution model, and then use it in a tree. )
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
We implementated support for rate variation

input: alpha OR insert your own custom rates
* simple example

### Site-Specific Transition Probabilities
* input: list of matrices
* input (optional): settng gtr.pi = [r1, r2, . . . rN] for single matrix W is exactly equivalent to inserting a list of matrices [r1 * W, r2 * W, . . . rN * W] 
* input (optional): gtr.mu
* simple example


### Changed Files
* significant changes, maybe a little bit about general purpose before us \
In treetime directory: \
treeanc.py \
gtr_site_specific.py \
seq_utils.py \
asrv.py (created) \
In test directory: \
run_test.py \
test_anc.py (created) \
test_treetime.py \

### Full Example
how to test-run the project on sample inputs to get the sample output
 
### Related Tools
__________ for generating tree toplogies, __________ for generating MSAs, __________ for infering alpha, ________ for infering protein structual information
 

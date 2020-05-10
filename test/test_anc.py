import numpy as np
from Bio import Phylo, AlignIO
from io import StringIO
from collections import defaultdict
from treetime import seq_utils
from treetime import TreeAnc, GTR
from treetime import ASRV
import random
import sys

try:
    from itertools import izip
except ImportError:  # python3.x
    izip = zip
    
def calc_mutation_rate(t):
    # estimate mutation rate
    node_order = []
    for node in t.tree.get_terminals():
        node_order.append(node)
    
    inconsistent_count = 0
    count = 0
    
    for i in range(len(node_order) - 1):
        for j in range(i + 1, len(node_order)):
            u = node_order[i]
            v = node_order[j]
            for k in range(len(u.sequence)):
                if u.sequence[k] != v.sequence[k]:
                    inconsistent_count += 1
            count += len(u.sequence)
    
    return inconsistent_count / count
    

if __name__ == '__main__':
    path = '../data/'
    #path = '../data/fp_seqs/'
    #path = '../data/treetime_data/h3n2_na/'
    
    
  
    # instantiate treetime
    #myTree = TreeAnc(gtr='jtt92', tree=path+'raxml_tree.nwk', aln=path+'h3n2_aa_aln.fasta', verbose=0)
    #myTree = TreeAnc(gtr='jtt92', tree=path+'phyml_tree.nwk', aln=path+'test_phyml.fasta', mu=5.0, verbose=0)
    
    #myTree = TreeAnc(gtr='t92', tree=path+'h3n2_na_500.nwk',
    #                  aln=path+'h3n2_na_500.fasta', verbose=0)
    #myTree.infer_ancestral_sequences(infer_gtr=False, marginal=True)
    #print(myTree.tree.sequence_LH.mean())

      
    # test ASRV
    ex_properties = 'EEEEEEEBEEEBEBEBEBEEEBEEEEBEBEEEEEEEEEEBEEEBEEEBEEEEEBEBBBEBBBEBBEEEEEEBEEEEEEBEEEBEEEEEEEEEEEEEEEEEEEBEBEBEEEBEBEEEEBEEEBEBEBEEEEEEEEEBEEEEEEEEEEEEEEEEEEEEBEBEBEEEBEBEEEEEEEBEBEEEEEBEEEEEEEEEBBBEEEBEEEEEEEEEEEEEEEEEBEBEEEEEE'
    solvent_property = list(ex_properties)
    myTree2 = TreeAnc(compress=True, gtr='wag', tree=path+'fp_tree.nwk', aln=path+'fp_extant_aa.fasta', verbose=0)
    myTree2.infer_ancestral_sequences(infer_gtr=False, marginal=True)
    
    test = myTree2.ancestral_likelihood()
    print(test.sum())
    
    print(myTree2.tree.sequence_LH.sum())
    
    # test ASRV + Mixture model
    #pseudo_solvent_accessibility = [np.random.choice(['B', 'E']) for i in range(469)]
    #myTree3= TreeAnc(gtr="EX", tree=path+'raxml_tree.nwk', aln=path+'h3n2_aa_aln.fasta',
    #                alpha=0.299, struct_propty=pseudo_solvent_accessibility, verbose=0)
    
    #myTree3.infer_ancestral_sequences(infer_gtr=False, marginal=True)
    #print(myTree3.tree.sequence_LH.mean())
    
    # test joint reconstruction
    #tree_joint = TreeAnc(gtr='jtt92', tree=path+'raxml_tree.nwk', aln=path+'h3n2_aa_aln.fasta', verbose=0)
    #tree_joint.infer_ancestral_sequences(infer_gtr=False, marginal=False)
    
    #tree_joint_asrv = TreeAnc(gtr='jtt92', tree=path+'raxml_tree.nwk', aln=path+'h3n2_aa_aln.fasta', alpha=0.299, verbose=0)
    #tree_joint_asrv.infer_ancestral_sequences(infer_gtr=False, marginal=False)
    
    #print(tree_joint.tree.sequence_LH.mean(), tree_joint_asrv.tree.sequence_LH.mean())
    

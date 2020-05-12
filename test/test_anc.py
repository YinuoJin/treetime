import numpy as np
from Bio import Phylo, AlignIO, SeqIO
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


def rename_ancestral_nodes(curr_tree):
    for node in curr_tree.tree.find_clades(order='postorder'):
        if len(node.clades) == 0:
            continue
        children_set = set([child.name for child in node.clades])
        if '01' in children_set and '02' in children_set:
            node.name = '37'
        elif '03' in children_set and '04' in children_set:
            node.name = '36'
        elif '05' in children_set and '06' in children_set:
            node.name = '35'
        elif '07' in children_set and '08' in children_set:
            node.name = '33'
        elif '09' in children_set and '10' in children_set:
            node.name = '31'
        elif '12' in children_set and '13' in children_set:
            node.name = '29'
        elif '14' in children_set and '15' in children_set:
            node.name = '27'
        elif '16' in children_set and '17' in children_set:
            node.name = '26'
        elif '18' in children_set and '19' in children_set:
            node.name = '25'
        elif '25' in children_set and '26' in children_set:
            node.name = '24'
        elif '24' in children_set and '27' in children_set:
            node.name = '23'
        elif '11' in children_set and '29' in children_set:
            node.name = '28'
        elif '23' in children_set and '28' in children_set:
            node.name = '22'
        elif '35' in children_set and '36' in children_set:
            node.name = '34'
        elif '33' in children_set and '34' in children_set:
            node.name = '32'
        elif '31' in children_set and '32' in children_set:
            node.name = '30'
        elif '22' in children_set and '30' in children_set:
            node.name = '21'
        elif '21' in children_set and '37' in children_set:
            node.name = '20'


def measure_accuracy(curr_tree, show_results=False):
    rename_ancestral_nodes(curr_tree)
    node_seq_dict = {}
    
    for node in curr_tree.tree.find_clades(order='preorder'):
        if node.name != 20:
            node_seq_dict[node.name] = node.sequence

    # import ground truth sequences
    true_node_seq_dict = {}
    record_dict = SeqIO.to_dict(SeqIO.parse(path + "fp_seq_aa.fasta", "fasta"))
    
    for name in node_seq_dict.keys():
        if name != '20':
            true_node_seq_dict[name] = np.array(record_dict[name].seq, dtype='U1')
    
    name_err = {}
    total_count = 0
    for name in true_node_seq_dict.keys():
        incorrect_count = ((node_seq_dict[name] == true_node_seq_dict[name]) == False).astype(int).sum()
        if int(name) > 20:
            name_err[name] = incorrect_count
            
        total_count += incorrect_count
        
    acc = 1 - total_count / (17 * len(curr_tree.tree.root.sequence))
    
    if show_results:
        print(total_count)
        print('acccuracy:', 1 - total_count / (17 * len(curr_tree.tree.root.sequence)))
    
    return total_count, acc, name_err


def sort_dict(dic):
    sorted_dict = []
    for key in sorted(dic.keys()):
        sorted_dict.append(dic[key])
    
    return sorted_dict


if __name__ == '__main__':
    path = '../data/'

    # test ASRV
    ex_properties = 'EEEEEEEBEEEBEBEBEBEEEBEEEEBEBEEEEEEEEEEBEEEBEEEBEEEEEBEBBBEBBBEBBEEEEEEBEEEEEEBEEEBEEEEEEEEEEEEEEEEEEEBEBEBEEEBEBEEEEBEEEBEBEBEEEEEEEEEBEEEEEEEEEEEEEEEEEEEEBEBEBEEEBEBEEEEEEEBEBEEEEEBEEEEEEEEEBBBEEEBEEEEEEEEEEEEEEEEEBEBEEEEEE'
    solvent_property = list(ex_properties)
    
    from treetime import ASRV
    
    asrv = ASRV(0.545)
    rates = asrv.calc_rates()
    myTree = TreeAnc(compress=False, gtr='wag', tree=path+'fp_tree.nwk', rates=rates, aln=path+'fp_extant_aa.fasta', verbose=0)
    myTree.infer_ancestral_sequences(infer_gtr=False, marginal=True)
    
    
    
    print(myTree.tree.sequence_LH)  # print LH
    print(myTree.tree.total_sequence_LH)  # print total LH

    # calculate accuracy against the ground truth
    _ = measure_accuracy(myTree, show_results=True)

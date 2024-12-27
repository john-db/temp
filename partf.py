from decimal import Decimal
import numpy as np
import treeswift as ts
import pandas as pd
from tree_scorer import calc_prob, pf_cond_on_one_tree, cell_lineage_tree_prob
import argparse
from collections import Counter

def main(args):
    trees = args.pathtrees
    path = args.patherror
    alpha=args.alpha
    beta=args.beta
    cells = sorted(args.cells.split('#'))
    mutation = args.mutation
    num_samples = args.num_samples

    df = pd.read_csv(path, sep="\t", index_col=[0]).sort_values(by=["cell_id_x_mut_id"])

    I_mtr = df.values
    t1 = I_mtr * (1 - beta) / (alpha + 1 - beta)
    t2 = (1 - I_mtr) * beta / (beta + 1 - alpha)
    P = t1 + t2
    P[I_mtr == 3] = 0.5

    
    my_cell = []
    for i in range(len(df.index)):
        if df.index[i] in cells:
            my_cell += [i]
    cond_c = np.zeros(P.shape[0], dtype=np.int8)
    cond_c[my_cell] = 1

    all_cells = list(df.index)
    clade = {cell : 0 for cell in all_cells}
    clade.update({c : 1 for c in cells})

    cond_m = np.where(df.columns == mutation)[0][0]

    idx_to_cells = df.index
    
    cells_to_idx = {idx_to_cells[i]:i for i in range(len(idx_to_cells))}

    idxp1_to_cells = {str(i + 1) : idx_to_cells[i] for i in range(len(idx_to_cells))}

    trees_list = []
    tree_probs = []
    cond_probs = []
    tree_counts = []

    with open(trees,"r") as file:
        for line in file:
            trees_list += [line]
    
    tuples = []
    trees_list_res = trees_list[0:num_samples]
    unique_trees = list(set(trees_list_res))
    counts = Counter(trees_list_res)
    for t in unique_trees:
        # subtrees = [np.zeros(len(cells_to_idx), dtype=int)]
        tree = ts.read_tree_newick(t)
        subtrees = []
        for node in tree.traverse_preorder():
            subtree = np.zeros(len(cells_to_idx), dtype=np.int8)
            # subtree = [0] * len(cells_to_idx)
            leaves = [cells_to_idx[leaf.get_label()] for leaf in node.traverse_leaves()]
            subtree[leaves] = 1
            # if subtree not in subtrees:
            subtrees += [subtree]
        prob = calc_prob(P, subtrees, 0)
        # prob = cell_lineage_tree_prob(P, subtrees)
        tree_probs += [prob]
        cond_probs += [pf_cond_on_one_tree(P, subtrees, cond_c, cond_m)]
        tree_counts += [counts[t]]
        tuples += [(prob, tree.newick())]
        
    numerator = Decimal(0)
    denominator = Decimal(0)
    for i in range(len(tree_probs)):
        numerator += Decimal(tree_probs[i] * cond_probs[i][0] / cond_probs[i][1]) * Decimal(tree_counts[i])
        denominator += Decimal(tree_probs[i]) * Decimal(tree_counts[i])
    
    output = path + "\t" + trees + "\t" + mutation + "\t" + ",".join(cells) + "\t" + str(alpha) + "\t" + str(beta) + "\t" + str(num_samples) + "\t" + str(Decimal(numerator)) + "\t" + str(Decimal(denominator)) + "\t" + str(Decimal(numerator) / Decimal(denominator))
    print(output)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='run.py')

    parser.add_argument("-pe", "--patherror", type=str,                                                        
                        help="input genotype matrix", required=True)
    parser.add_argument("-pt", "--pathtrees", type=str,                                                        
                        help="sample of trees (SCITE labelled)", required=True)
    parser.add_argument("-n", "--num_samples", type=int,                                                        
                        help="num_samples", required=True)
    parser.add_argument("-C", "--cells", type=str,                                                        
                        help="cells (# sep)", required=True)
    parser.add_argument("-m", "--mutation", type=str,                                                        
                        help="mutation", required=True)
    parser.add_argument("-a", "--alpha", type=float,                                                        
                        help="alpha (fp)", required=True)
    parser.add_argument("-b", "--beta", type=float,                                                        
                        help="beta (fn)", required=True)
    main(parser.parse_args())
import numpy as np
import pandas as pd
import treeswift as ts
from decimal import Decimal
import numpy.linalg as la
import itertools
from cf_mat_to_newick import cf_to_newick

def comb(max, k):
    return [comb for comb in itertools.combinations(list(range(max)), k)]

def calc_prob(P, subtrees, order):
    nontrivial = list(subtrees[24:-1])
    sub_list = list(subtrees)

    logP = np.log2(P)
    logPC = np.log2(1 - P)
    p0 = 2 ** Decimal(prob_mat_mul_calc(logP, logPC, subtrees))
    correction = Decimal(0)
    for size in range(1, order + 1):
        temp = Decimal(0)
        for c in comb(22, size):
            restricted_subtrees = nontrivial[0:c[0]]
            for i in range(len(c)):
                if i < len(c) - 1:
                    restricted_subtrees += nontrivial[c[i] + 1 : c[i + 1]]
            restricted_subtrees += nontrivial[c[-1] + 1:]
            if restricted_subtrees != []:
                restricted_subtrees += subtrees[0:23]
                restricted_subtrees += [subtrees[-1]]
                temp += 2 ** Decimal(prob_mat_mul_calc(logP, logPC, np.array(restricted_subtrees)))
        if temp != 1:
            correction += (Decimal(-1) ** (size)) * temp
    return p0 + correction
            


def pf_cond_on_one_tree(P, subtrees, cond_c, cond_m):
    r"""
    Prob_{A\sim P}[\subtree(c, R, A)\cap A\in G| A\in T] in O(n^2).

    :param P:
    :param subtrees: cell lineage tree  n x (2n+1)
    :param cond_c: set of cells
    :param cond_m: one mutation
    :return: conditioned on the given tree what is probability of the given partition
        based on P numerator, denominator are returned separately here
    """

    denominator = Decimal(0)
    numerator = Decimal(0)
    col = P[:, cond_m]
    for v in subtrees:
        prob = Decimal(np.prod(col * v + (1 - col) * (1 - v)))
        denominator += prob
        if np.array_equal(v, cond_c):
            numerator = prob
    return numerator, denominator

def prob_mat_mul_calc(logP, logPC, subtrees):
    st = np.array(subtrees)

    res = np.exp2(np.matmul(st, logP) + np.matmul(1 - st, logPC))
    res = np.matmul(np.ones(res.shape[0]), res)
    res = np.matmul(np.ones(res.shape[0]), np.log2(res))

    return res

def prob_mat_vec_calc(P, subtrees):
    logP = np.log2(P)
    logPC = np.log2(1 - P)
    st = np.array(subtrees)
    # res = np.exp2(np.matmul(st, lp))

    return_value = Decimal(1)
    for j in range(P.shape[1]):
        # temp = Decimal(0)
        colP = logP[:, j]
        colPC = logPC[:, j]
        res = np.exp2((np.matmul(st, colP) + np.matmul(1 - st, colPC)))
        temp = Decimal(np.matmul(np.ones(res.shape[0]), res))
        return_value *= temp
    return return_value

def cell_lineage_tree_prob(P, subtrees):
    r"""
    Calculate Prob_{A\sim P}[A\in T] in O(n m^2).

    :param P:
    :param subtrees: cell lineage tree
    :return: Probability of this tree in the original distribution( based on P)
    """

    # subtrees.sort(key=sum)
    # for tree in subtrees:
    #     print(tree)

    return_value = Decimal(1)
    for j in range(P.shape[1]):
        temp = Decimal(0)
        col = P[:, j]
        for v in subtrees:
            temp += Decimal(np.prod(col * v + (1 - col) * (1 - v)))
        return_value *= temp
    return return_value

def main():
    restrict = 24

    path = "/Users/john/Desktop/matrices/Sep9Matrices/D-C19_and_C2C5_corrected.tsv"
    path = "/Users/john/Desktop/matrices/Sep9Matrices/D-C19_corrected.tsv"
    df = pd.read_csv(path, sep="\t", index_col=[0]).sort_values(by=["cell_id_x_mut_id"])
    df = df.head(restrict)

    # path_cf = "/Users/john/Desktop/e_data/E-jun_13_2024-a_1e-8-b_0.075.tsv"
    path_cf = "/Users/john/Desktop/matrices/Sep9Matrices/E-C19_and_C2C5_corrected-fp_0.0001-fn_0.075.tsv"
    path_cf = "/Users/john/Desktop/matrices/Sep9Matrices/E-C19_corrected-fp_0.001-fn_0.075.tsv"
    cf = pd.read_csv(path_cf, sep="\t", index_col=[0]).sort_values(by=["cell_id_x_mut_id"])
    cf = cf.head(restrict)


    newick = cf_to_newick(cf)
    # newick = "(C13,(C9,(((C10,C12),((C23,(C24,(C6,(C19,C21)))),(C14,C3))),(((C1,(C22,C4)),(C17,(C2,C5))),((C18,(C15,(C8,(C20,C7)))),(C11,C16))))))"    


    alpha=0.0001
    beta=0.075



    idx_to_cells = df.index
    cells_to_idx = {idx_to_cells[i]:i for i in range(len(idx_to_cells))}


    I_mtr = df.values
    t1 = I_mtr * (1 - beta) / (alpha + 1 - beta)
    t2 = (1 - I_mtr) * beta / (beta + 1 - alpha)
    P = t1 + t2
    P[I_mtr == 3] = 0.5

    tree = ts.read_tree_newick(newick)
    # subtrees = [np.zeros(len(cells_to_idx), dtype=int)]
    subtrees = []
    for node in tree.traverse_preorder():
        subtree = np.zeros(len(cells_to_idx), dtype=np.int8)
        # subtree = [0] * len(cells_to_idx)
        leaves = [cells_to_idx[leaf.get_label()] for leaf in node.traverse_leaves()]
        subtree[leaves] = 1
        # if subtree not in subtrees:
        subtrees += [subtree]

    subtrees = [np.array(x) for x in {(tuple(e)) for e in subtrees}]
    subtrees.sort(key=sum)

    print(cell_lineage_tree_prob(P, subtrees))
    for order in range(0, 24):
        p = calc_prob(P, subtrees, order)
        print("order " + str(order) + " " + str(p))

    

if __name__ == "__main__":
    main()
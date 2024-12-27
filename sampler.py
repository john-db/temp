import pandas as pd
import argparse
import _partition_function
from _partition_function import partition_function
from decimal import Decimal
import os
import pathlib
import treeswift as ts
import numpy as np

from _clt_sampler import draw_sample_clt
from cf_mat_to_newick import cf_to_newick

def main(args):
    path = args.input_matrix
    df = pd.read_table(path, index_col=0)
    
    alpha=args.alpha
    beta=args.beta

    divide = True
    num_samples = args.num_samples
    delta = 0.8
    eps = 10.0
    coef = 10
    seed = args.seed

    output = "number of samples = " + str(num_samples) + "\n"
    cells = sorted(args.cells.split(","))

    muts = [args.mutation]
    names_to_cells = list(df.index)

    external_folder = str(pathlib.Path(__file__).parent.resolve()) + "/external/"

    os.system('rm ' + external_folder + './matrix_for_scite*')
    os.system('rm ' + external_folder + 'PARTF_ACCEPTED_NEWICKS.temp')
    pd.DataFrame(df.transpose().values).to_csv(external_folder + "matrix_for_scite.temp", header=False, index=False, sep=" ")
    cmd = external_folder + 'scite -i ' + external_folder + 'matrix_for_scite.temp -n ' + str(len(df.columns)) + ' -m ' + str(len(df)) + ' -r 1 -l ' + str(args.mcmc_length) + ' -fd ' + str(alpha) + ' -ad ' + str(beta) + ' -transpose -newicks ' + external_folder + 'PARTF_ACCEPTED_NEWICKS.temp'
    print(cmd)
    os.system(cmd)

    mcmc_trees = []
    idxp1_to_cells = {str(i + 1) : df.index[i] for i in range(len(df.index))}
    with open(external_folder + 'PARTF_ACCEPTED_NEWICKS.temp', 'r') as file:
        for line in file:
            newick = line
            tree = ts.read_tree_newick(newick + ";")
            tree.rename_nodes({label : None for label in tree.labels(leaves=False, internal=True)})
            tree.rename_nodes(idxp1_to_cells)
            mcmc_trees += [tree.newick()]


    I_mtr = df.values
    t1 = I_mtr * (1 - beta) / (alpha + 1 - beta)
    t2 = (1 - I_mtr) * beta / (beta + 1 - alpha)
    P = t1 + t2
    P[I_mtr == 3] = 0.5
    P = P.astype(np.float64)

    rng = np.random.default_rng(seed=args.seed)
    sampled_trees = []
    probs = []
    for i in range(args.num_samples):
        rand = rng.random()
        tree = None
        if rand < args.split:
            _, subtrees, prior_prob = draw_sample_clt(P, False, eps=eps, delta=delta, divide=divide, coef=coef, names_to_cells=names_to_cells, clade=cells, rng=rng)
            arr = np.array(subtrees).T
            # mat = np.array([s for s in subtrees if sum(s) > 1]).T
            newick = cf_to_newick(pd.DataFrame(arr, index=df.index))
            tree = newick + ";"
        else:
            choice = rng.integers(low=0, high=len(mcmc_trees))
            tree = mcmc_trees[choice]
        sampled_trees += [tree]
    # with open(args.output, 'w') as f:
    #     for tree in sampled_trees:
    #         f.write(f"{tree}\n")
    os.system('rm ' + external_folder + './matrix_for_scite*')
    os.system('rm ' + external_folder + 'PARTF_ACCEPTED_NEWICKS.temp')
    for tree in sampled_trees:
        print(tree)
        

        

    # pf = partition_function(df_input=df, alpha=alpha, beta=beta, n_samples=num_samples, n_batches=1, muts=muts, cells=cells, names_to_cells=names_to_cells,eps = eps, delta=delta, divide=divide, coef=coef, mcmc_trees=mcmc_trees, split_prob=args.split, my_seed=seed)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='run.py')

    parser.add_argument("-i", "--input_matrix", type=str,                                                        
                        help="Path to input genotype matrix where rows correspond to cells/sublines and columns correspond to mutations. See repo examples for formatting.", required=True)
    parser.add_argument("-n", "--num_samples", type=int,                                                        
                        help="Number of trees to be sampled", required=True)
    parser.add_argument("-c", "--cells", type=str,                                                        
                        help="List of cells (comma separated)", required=True)
    parser.add_argument("-m", "--mutation", type=str,                                                        
                        help="Name of the mutation (column) in the matrix", required=True)
    parser.add_argument("-fp", "--alpha", type=float,                                                        
                        help="False-positive rate (alpha in the paper)", required=True)
    parser.add_argument("-fn", "--beta", type=float,                                                        
                        help="False-negative rate (beta in the paper)", required=True)
    parser.add_argument("-sp", "--split", type=float,                                                        
                        help="Split probability (probability of picking a tree by independent sampling).", required=False, default=0.5)
    parser.add_argument("-l", "--mcmc_length", type=int,                                                        
                        help="How many MCMC trees to collect.", required=False, default=10000)
    parser.add_argument("-s", "--seed", type=int,                                                        
                        help="random seed", required=False, default=0)
    
    # -i "/Users/john/Desktop/journal_partf/Partition-Function/partition_function/example/data.tsv" -n 2 -c "cell1,cell2" -m "mut1" -fp 0.01 -fn 0.1 -sp 0.5 -l 10 -s 0
    main(parser.parse_args())

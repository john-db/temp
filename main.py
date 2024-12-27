import pandas as pd
import argparse
import _partition_function
from _partition_function import partition_function
from decimal import Decimal
import os
import pathlib
import treeswift as ts

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

    pd.DataFrame(df.transpose().values).to_csv(external_folder + "matrix_for_scite.temp", header=False, index=False, sep=" ")
    
    os.system(external_folder + 'scite -i ' + external_folder + 'matrix_for_scite.temp -n ' + str(len(df.columns)) + ' -m ' + str(len(df)) + ' -r 1 -l ' + str(args.mcmc_length) + ' -fd ' + str(alpha) + ' -ad ' + str(beta) + '-transpose -newicks ' + external_folder + 'PARTF_ACCEPTED_NEWICKS.temp')

    mcmc_trees = []
    idxp1_to_cells = {str(i + 1) : df.index[i] for i in range(len(df.index))}
    with open(external_folder + 'PARTF_ACCEPTED_NEWICKS.temp', 'r') as file:
        for line in file:
            newick = line
            tree = ts.read_tree_newick(newick + ";")
            tree.rename_nodes(idxp1_to_cells)
            tree.rename_nodes({label : None for label in tree.labels(leaves=False, internal=True)})
            mcmc_trees += [tree.newick()]

    pf = partition_function(df_input=df, alpha=alpha, beta=beta, n_samples=num_samples, n_batches=1, muts=muts, cells=cells, names_to_cells=names_to_cells,eps = eps, delta=delta, divide=divide, coef=coef, mcmc_trees=mcmc_trees, split_prob=args.split, my_seed=seed)
    
    os.system('rm ' + external_folder + './matrix_for_scite*')
    os.system('rm ' + external_folder + 'PARTF_ACCEPTED_NEWICKS.temp')


    num = pf[0].iloc[0][0]
    denom = pf[0].iloc[0][1]


    output = "Input: " + args.input_matrix + "\n"
    output += "False-positive rate: " + str(args.alpha) + "\n"
    output += "False-negative rate: " + str(args.beta) + "\n"
    output += "Clade: " + ",".join(cells) + "\n"
    output += "Mutation: "  + args.mutation + "\n"
    output += "Number of samples: " + str(args.num_samples) + "\n"
    output += "RNG seed: " + str(seed) + "\n"
    output += "numerator_est: " + str(num) + "\n"
    output += "denominator_est: " + str(denom) + "\n"
    output += "Partition function estimate: " + str(float(Decimal(num / denom)))
    print(output)



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
    parser.add_argument("-sp", "--split", type=float,                                                        
                        help="Split probability (probability of picking a tree by independent sampling).", required=False, default=0.5)
    parser.add_argument("-l", "--mcmc_length", type=int,                                                        
                        help="How many MCMC trees to collect.", required=False, default=10000)
    parser.add_argument("-fn", "--beta", type=float,                                                        
                        help="False-negative rate (beta in the paper)", required=True)
    parser.add_argument("-s", "--seed", type=int,                                                        
                        help="random seed", required=False, default=0)
    
    
    main(parser.parse_args())

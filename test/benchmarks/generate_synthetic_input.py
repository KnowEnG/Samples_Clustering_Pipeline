"""
Created on Mon Sep 8 10:32:25 2016
@author: Xi
@author: The Gene Sets Characterization dev team
"""
import itertools
import numpy as np
import pandas as pd

def generate_synthetic_data(row_num, col_num, k):
    """This is to generate synthetic data used in the
    benchmark test.

    Args:
        row_num: number of rows in user spreadsheet
        col_num: number of cols in user spreadsheet
        k: k clusters
    """
    g_name = ['g'+str(i) for i in np.arange(row_num)]
    u_name = ['user'+str(i) for i in np.arange(col_num)]
    g = len(g_name)/k
    u = len(u_name)/k
    spreadsheet_df = pd.DataFrame(columns=u_name, index=g_name).fillna(0)

    for i in range(k):
        spreadsheet_df.values[i*g:(i+1)*g, i*u:(i+1)*u] = 1
    spreadsheet_df.to_csv("synthetic_user_spreadsheet.df", header=True, index=True, sep='\t')

    val = []
    gene_names = spreadsheet_df.index.values
    for i in range(k):
        combine_two = itertools.combinations(gene_names[g*i: g*(i+1)], 2)
        for i in combine_two:
            val.append(list(i)+[1, "test"])
    edge_file = pd.DataFrame(val, columns=['g1', 'g2', 'wt', 'type'], index=None)
    edge_file.to_csv("synthetic_gg_network.edge", header=False, index=False, sep='\t')

generate_synthetic_data(420, 120, 3)

import numpy as np
import pandas as pd
import itertools

def generate_benchmark_data(row_num, col_num, k):
    g_name = ['g'+str(i) for i in np.arange(row_num)]
    u_name = ['user'+str(i) for i in np.arange(col_num)]
    g = len(g_name)/k
    u = len(u_name)/k
    spreadsheet_df = pd.DataFrame(columns=u_name, index=g_name).fillna(0)
    
    for i in range(k):
        spreadsheet_df.values[i*g:(i+1)*g, i*u:(i+1)*u] = 1
    spreadsheet_df.to_csv("user_spreadsheet_benchmark.txt", header=True, index=True, sep='\t')

    val = []
    li = spreadsheet_df.index.values
    for i in range(k):
        a = itertools.combinations(li[g*i: g*(i+1)], 2)
        for i in a:
            val.append(list(i)+[1,"test"])
    edge_file = pd.DataFrame(val, columns=['g1', 'g2', 'wt', 'type'], index=None)
    edge_file.to_csv("gg_network_benchmark.txt", header=False, index=False, sep='\t')
    
generate_benchmark_data(420, 120, 3)

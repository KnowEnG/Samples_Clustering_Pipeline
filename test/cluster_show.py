"""
display clustering
test tool expecing 2 column .tsv with cluster assignments in second column

lanier4@illinois.edu
"""
import sys
import argparse
import pandas


def main(args):
	parser = argparse.ArgumentParser()
	parser.add_argument('-file_name', type=str)
	args = parser.parse_args()
	f_name = args.file_name
	clusters_df = pandas.read_csv(f_name, sep='\t', header=None, index_col=0)
	for k in range(0, clusters_df[1].max()+1):
		print(k, (clusters_df[1] == k).sum())

if __name__ == "__main__":
	main(sys.argv)
#! /usr/bin/env python2.7

from optparse import OptionParser
import os, sys
import PCA_Utils as util

"""
subsample_dataset.py
Author: Alex Sevy <alex.sevy@gmail.com>

This script is used to generate subsampled replicates of IGHV-IGHJ gene 
frequency data from antibody repertoire sequencing. It is the first step
of the repertoire fingerprinting protocol described in Sevy, Soto et al. 2018
"""



if __name__ == '__main__':
	usage = 'Usage: python subsample_datasets.py --reps=10 --depth=100000 --n_jobs=8 hip1.tsv hip2.tsv'
	parser=OptionParser(usage)

	parser.add_option('--reps',dest='reps', help='How many simulated datasets to make? Default=10', default=10)
	parser.add_option('--depth',dest='depth', help='Subsample dataset to what depth? Default=10^5', default=10**5)
	parser.add_option('--n_jobs',dest='n_jobs', help='How many processes should be run? Default=1', default=1)

	(options,args)= parser.parse_args()

	## If no arguments provided, print usage and exit gracefully
	if len(args) == 0:
		print usage
		exit()

	## All output data will be stored in a folder called datasets/
	if not os.path.exists('datasets'):
		os.mkdir('datasets')

	for dataset in args:
		print 'working on %s'%dataset

		## read_table eliminates all but the 306 prefined features, and adds
		## back 0 values for V-J pairs that aren't observed
		df = util.read_table(dataset)

		## run subsampling
		subsampled_dfs = util.subsample( df, reps=int(options.reps), size=int(options.depth), n_jobs=int(options.n_jobs) )
		
		## output the subsampled dataframes
		for ii, new_df in enumerate( subsampled_dfs ):
			print 'writing number %d'%(ii+1)
			out_name = 'datasets/%s_%d.tsv'%(dataset.split('.')[0],ii+1)
			new_df.to_csv( out_name, index=True, sep='\t' ) 
#! /usr/bin/env python2.7

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import PCA_Utils as util
import collections
from sklearn.decomposition import PCA
import seaborn as sns
from optparse import OptionParser
import os

script_path = os.path.dirname(os.path.realpath(__file__))
style_file = os.path.join( script_path, 'style', 'prism_color.mplstyle' )
plt.style.use(style_file)

"""
analyze_pca.py
Author: Alex Sevy <alex.sevy@gmail.com>

This script is used to run principal component analysis on  IGHV-IGHJ gene 
frequency data from antibody repertoire sequencing. It is the second step
of the repertoire fingerprinting protocol described in Sevy, Soto et al. 2018,
and is meant to be preceded by subsample_datasets.py
"""


if __name__ == '__main__':
	usage = 'Usage: python analyze_pca.py hd1 hd2 hd3 ....'
	parser=OptionParser(usage)

	parser.add_option('--output',dest='output', help='Output prefix for PDFs. Default=pca', default='pca')
	parser.add_option('--ncomp',dest='ncomp', help='How many components for PCA. Default=8', default=8)
	parser.add_option('--reps',dest='reps', help='How many simulated datasets were made? Default=10', default=10)
	parser.add_option('--path',dest='path', help='Path where files are located. Default=datasets/', default='datasets')

	(options,args)= parser.parse_args()

	if len(args) == 0:
		print usage
		exit()


	data_dict = {}
	norm_dict = {}

	for key in args:

		# go through each subsampled replicates one by one, load data and normalize
		# by Z score
		for ii in range(1,int(options.reps)+1):
			data_dict['%s_%d'%(key,ii)] = util.read_table( '%s/%s_%d.tsv'%(options.path,key,ii) )
			norm_dict['%s_%d'%(key,ii)] = util.normalize_zscore( data_dict['%s_%d'%(key,ii)] )

			## Sort rows so they are in numeric not alphabetical order
			df = norm_dict[ '%s_%d'%(key,ii) ]
			rows = list(df.index)
			new_rows = map(lambda a: [int(aa) for aa in a[4:].split('-')], rows)
			new_rows.sort()
			new_row_labels = map(lambda a: '-'.join( ['IGHV%d'%a[0]]+[str(i) for i in a[1:]]), new_rows )

			df = df.reindex(new_row_labels)
			norm_dict[ '%s_%d'%(key,ii) ] = df
		
	# set which data points are to be used to train PCA (trainset) and which
	# should only have PCA applied to them (testset). Defaults to training
	# and testing on all data points
	trainset = args 
	testset = args

	labels = [item.split('/')[-1].split('.')[0] for item in args]

	# flatten the tables from a dataframe to a 1D array
	Xtrain = [norm_dict['%s_%d'%(key,ii)].values.flatten() for key in trainset for ii in range(1,int(options.reps)+1)]
	Xtest = [norm_dict['%s_%d'%(key,ii)].values.flatten() for key in testset for ii in range(1,int(options.reps)+1)]

	# train PCA
	ncomp = int(options.ncomp)
	pca = PCA(ncomp)
	pca.fit(Xtrain)

	# Make scree plot with explained variance for each component
	plt.figure()
	bars = plt.bar( range(1,ncomp+1), pca.explained_variance_ratio_, align='center' )
	ax = plt.gca()
	for rect in bars:
		height = rect.get_height()
		ax.text(rect.get_x() + rect.get_width()/2., 1.0*height,
			'%0.2f' % float(height),
			ha='center', va='bottom')

	plt.ylabel('Explained variance ratio')
	plt.xlabel( 'Principal component')
	plt.savefig( options.output+'_explained_variance.pdf' )
	plt.close()

	# Transform data and plot in 2D
	plt.figure()
	ax = plt.gca()

	for key in testset:
		tf = pca.transform(
			[norm_dict['%s_%d'%(key,ii)].values.flatten() for ii in range(1,int(options.reps)+1)]
			)

		x = tf[:,0]
		y = tf[:,1]

		plt.plot( x, y,'o', markersize=6 )

		ax.annotate(key, xy=(max(x),max(y)))

	plt.xlabel('PC1 (%0.0f%%)'%(pca.explained_variance_ratio_[0]*100))
	plt.ylabel('PC2 (%0.0f%%)'%(pca.explained_variance_ratio_[1]*100))
	plt.savefig(options.output+'_2D_projection.pdf')
	plt.close()


	cols = util.get_jgenes()
	rows = util.get_vgenes()

	## take a 1D array and turn it into a dataframe
	pc1_df = util.array_to_df( pca.components_[0], cols, rows)
	pc2_df = util.array_to_df( pca.components_[1], cols, rows)

	plt.figure(figsize=(5,10))
	sns.heatmap( pc1_df,cmap='RdBu_r',center=0,square=True, yticklabels=rows, xticklabels=cols)
	plt.yticks(rotation=0)
	plt.xticks(rotation=90)
	plt.title('PC1')
	plt.savefig(options.output+'_PC1_heatmap.pdf')
	plt.close()

	plt.figure(figsize=(5,10))
	sns.heatmap( pc2_df,cmap='RdBu_r',center=0,square=True, yticklabels=rows, xticklabels=cols)
	plt.yticks(rotation=0)
	plt.xticks(rotation=90)
	plt.title('PC2')
	plt.savefig(options.output+'_PC2_heatmap.pdf')
	plt.close()

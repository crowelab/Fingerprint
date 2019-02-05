#! /usr/bin/env python2.7

import pandas as pd
import numpy as np
import math
import multiprocessing as mp
from functools import partial
from tqdm import tqdm

def subsample( df, reps=10, size=1e5, n_jobs=1 ):
	pool = mp.Pool( processes=n_jobs )

	make_subsample_partial = partial(make_subsample, df=df, size=size)

	samples = list(tqdm(
		pool.imap( make_subsample_partial, range(int(reps)) ),
			total=int(reps) ))

	assert len(samples)==reps
	if reps > 1:
		assert not (samples[0]==samples[1]).all().all()
	return samples

def get_keys( df ):
	rows = list(df.index)
	cols = list(df)
	return [a+':'+b for a in rows for b in cols]

def make_subsample(i, df, size):
	np.random.seed()
	labels = np.array(get_keys( df ))
	array = df.values.flatten()
	frequencies = array/float(array.sum())
	sampled_df = pd.DataFrame(0, columns=list(df),index=list(df.index),dtype=np.int)
	for item in np.random.choice(labels, int(size), p=frequencies):
		v,j = item.split(':')
		sampled_df[j][v] += 1
	return sampled_df

def get_vgenes():
	return map(lambda a: 'IGHV'+a, ['1-2','1-3','1-8','1-18','1-24','1-45','1-46','1-58','1-69','1-69-2','2-5','2-26','2-70','3-7','3-9','3-11','3-13','3-15','3-20','3-21','3-23','3-30','3-30-3','3-33','3-43','3-47','3-48','3-49','3-52','3-53','3-64','3-66','3-71','3-72','3-73','3-74','4-4','4-28','4-30-2','4-30-4','4-31','4-34','4-38-2','4-39','4-55','4-59','4-61','5-10-1','5-51','6-1','7-4-1'])

def get_jgenes():
	return map(lambda a: 'IGHJ'+a, ['1','2','3','4','5','6'])

def read_table( filename ):
	vgenes = get_vgenes()
	jgenes = get_jgenes()

	df = pd.read_table(filename, index_col=0, delim_whitespace=True)
	drop_genes = [col for col in list(df.index) if col not in vgenes]
	df = df.drop(drop_genes)

	# add back missing genes
	add_jgenes = [col for col in jgenes if col not in list(df)]
	add_vgenes = [col for col in vgenes if col not in list(df.index)]

	for gene in add_jgenes:
		df[gene] = [0]*len(list(df.index))

	for gene in add_vgenes:
		s = pd.Series([0]*6, index=df.columns)
		s.name = gene
		df = df.append( s )
	
	df = df.sort_index(axis=1)
	df = df.sort_index(axis=0)
	return df

def normalize_zscore( df, pseudocount=0.01 ):
	df = df.replace(0,pseudocount).apply(np.log10)
	norm = (df-df.values.mean())/df.values.std()
	return norm.apply(lambda a: 10**a)

def array_to_df( array, cols, rows ):
	pc1_df = pd.DataFrame(0, index=rows, columns=cols)
	for ii, component in enumerate(array):
		pc1_df.iloc[int(math.floor(ii/len(cols))),ii%len(cols)] = component
	return pc1_df

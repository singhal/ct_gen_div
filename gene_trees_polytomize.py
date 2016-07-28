import glob
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import os
import re

dir = '/Volumes/heloderma4/sonal/eco_IBD_oz/data/species_tree/gene_trees_Ct/'
files = glob.glob('%s*bestTree.tre' % dir)

ape = importr('ape')
for file in files:
	newfile = file.replace('bestTree', 'bestTree_poly')
	a = ape.read_tree(file)
	a = ape.di2multi(a, tol=5e-6)
	ape.write_tree(a, file=newfile)

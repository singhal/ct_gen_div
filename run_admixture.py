import re 
import glob
import pandas as pd
import os
import subprocess

c_file1 = '/Volumes/heloderma4/sonal/eco_IBD_oz/data/clustering/Ctenotus_Lerista.clustering.revised.csv'
d = pd.read_csv(c_file1)
d = d[d.GMYC_RAxML2.notnull()]

clusters = dict([(name, group.sample.tolist()) for name, group in d.groupby('GMYC_RAxML2')])


for cl, inds in clusters.items():
	if len(inds) > 2:
		maxK = 5
		if len(inds) < maxK:
			maxK = int(len(inds) / 2.0)
			if maxK < 2:
				maxK = 2
		for i in range(1, maxK + 1):
			out = '%s.%s.Q' % (cl, i)
			if not os.path.isfile(out):
				print '/Volumes/heloderma4/sonal/bin/admixture_macosx-1.23/admixture --cv %s.geno %s -j4 > %s.K%s.log' % (cl, i, cl, i)

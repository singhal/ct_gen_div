import re 
import glob
import pandas as pd
import os
import subprocess

c_file1 = '/Volumes/heloderma4/sonal/eco_IBD_oz/data/clustering/Ctenotus_Lerista.clustering.revised.csv'
d = pd.read_csv(c_file1)
d = d[d.GMYC_RAxML2.notnull()]

dir = '/Volumes/heloderma4/sonal/eco_IBD_oz/data/admixture/'

clusters = dict([(name, group.sample.tolist()) for name, group in d.groupby('GMYC_RAxML2')])

print 'cluster,inds,K'
for cl, inds in clusters.items():
        if len(inds) > 2:
		logs = glob.glob('%s%s.K*.log' % (dir, cl))
		cv = {}
		for log in logs:
			k = int(re.search('K(\d+)', log).group(1))
			f = open(log, 'r')
			for l in f:
				if re.search('CV', l):
					val = re.search('([0-9|\.]+)$', l.rstrip()).group(1)
					cv[k] = float(val)
		keys = sorted(cv.keys())
		winner = keys[0]
		cur_cv = None
		for key in keys:
			if cur_cv:
				if cv[key] < cur_cv:
					winner = key
					cur_cv = cv[key]
				else:
					break
			else:
				cur_cv = cv[key]
		print '%s,%s,%s' % (cl, len(inds), winner)

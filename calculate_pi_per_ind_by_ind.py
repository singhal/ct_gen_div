import re
import subprocess
import pandas as pd
import os
import numpy as np
import argparse
import gzip

parser = argparse.ArgumentParser(description="Calculate het and Fst.")
parser.add_argument('--ind', help="Cluster for which to run.")
parser.add_argument('--cov', help="coverage limit.")
args = parser.parse_args()
ind = args.ind
cov = args.cov
cov = int(args.cov)

c_file = '/Volumes/heloderma4/sonal/eco_IBD_oz/data/clustering/Ctenotus_Lerista.clustering.revised2.csv'
vcf_dir = '/Volumes/heloderma4/sonal/eco_IBD_oz/snp_calling/ind_variants/'
out_dir = '/Volumes/heloderma4/sonal/eco_IBD_oz/data/diversity/'


def get_clusters(c_file):
        d = pd.read_csv(c_file)
        d = d[d.GMYC_RAxML2.notnull()]

        inds = dict(zip(d['sample'], d.GMYC_RAxML2))

        return inds


def get_diversity(cl, ind, vcf_dir, out_dir):
	file = '%s%s.coverage_filtered_%s.vcf.gz' % (vcf_dir, ind, cov)

	# initialize counter
	hets = {'diff': 0, 'denom': 0}

	f = gzip.open(file, 'r')
	for l in f:
		if not re.search('#', l) and not re.search('INDEL', l):
                        d = re.split('\s+', l.rstrip())
                        # don't mess with multiallelics
                        if len(re.split(',', d[4])) == 1:
				geno = d[9]
				geno = re.search('^(\S\/\S)', geno).group(1)

				if geno in ['0/0', '0/1', '1/1']:
					hets['denom'] += 1
					if geno == '0/1':
						hets['diff'] += 1
	f.close()

	out_file = '%sall_clusters.individual_pi_coverage%s.csv' % (out_dir, cov)
	o = open(out_file, 'a')
	if hets['denom'] > 0:
		pi = hets['diff']  / float(hets['denom'])
	else:
		pi = np.nan
	o.write('%s,%s,%s,%.6f\n' % (cl, ind, hets['denom'], pi))
	o.close()


def initiate_file(out_dir):
	out_file = '%sall_clusters.individual_pi_coverage%s.csv' % (out_dir, cov)
	if not os.path.isfile(out_file):
		o = open(out_file, 'w')
		o.write('cluster,ind,denom,pi\n')
		o.close()


clusters = get_clusters(c_file)
initiate_file(out_dir)
get_diversity(clusters[ind], ind, vcf_dir, out_dir)

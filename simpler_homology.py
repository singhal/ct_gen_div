import re
import pandas as pd
import subprocess
import os
import argparse

parser = argparse.ArgumentParser(description="Run simple homology search.")
parser.add_argument('--cl', help="Cluster for which to run homology search.")
args = parser.parse_args()
cl = args.cl

# clustering file
c_file = '/Volumes/heloderma4/sonal/eco_IBD_oz/data/clustering/Ctenotus_Lerista.clustering.revised.csv'
# directory
dir = '/Volumes/heloderma4/sonal/eco_IBD_oz/snp_calling/rainbow_assemblies/'

WCLUST = 0.97


def get_clusters(c_file, cl):
        '''
        get the inds in a given cluster
        '''
        d = pd.read_csv(c_file)
        
        inds = d[d.GMYC_RAxML2 == cl].sample.tolist()
       
        return inds


def get_homolog(cluster, inds, dir):
	subdir = dir + 'indiv/'
	tmp1 = '%s%s.tmp1.fa' % (subdir, cluster)
	tmp2 = '%s%s.tmp2.fa' % (subdir, cluster)
	tmp3 = '%s%s.tmp3.fa' % (subdir, cluster)
	final = '%s%s.fa' % (dir, cluster)

	inds = ['%s%s.fa' % (subdir, ind) for ind in inds]
	new_inds = []

	for ind in inds:
		indtmp1 = ind + '1'
		indtmp2 = ind + '2'
		subprocess.call("vsearch --derep_fulllength %s --output %s --fasta_width 0" % (ind, indtmp1), shell=True)
        	subprocess.call("vsearch --cluster_smallmem %s --centroids %s --id %s --usersort --fasta_width 0" % (indtmp1, indtmp2, WCLUST), shell=True)
		new_inds.append(indtmp2)
		os.remove(indtmp1)

	subprocess.call("cat %s > %s" % (' '.join(new_inds), tmp1), shell=True)
	call = [os.remove(x) for x in new_inds]
	subprocess.call("vsearch --derep_fulllength %s --output %s --sizeout --fasta_width 0 --threads 4" % (tmp1, tmp2), shell=True)
	subprocess.call("vsearch --cluster_smallmem %s --centroids %s --sizein --sizeout --id %s --usersort --fasta_width 0 --threads 4" % (tmp2, tmp3, WCLUST), shell=True)

	# the number of times a contig should be sampled
	ninds = len(new_inds)
	if ninds < 3:
		cl_depth = 1
	elif ninds in [3, 4, 5]:
		cl_depth = 2
	else:
		cl_depth = round(ninds / 3.0)


	f = open(tmp3, 'r')
	o = open(final, 'w')
	ix = 0
	for l in f:
	        if re.search('>', l):
			seq = f.next().rstrip()
                	size = int(re.search('size=(\d+)', l).group(1))
                	if size >= cl_depth:
				o.write('>%s_%s\n%s\n'% (cluster, ix, seq))
				ix += 1
	f.close()
	o.close()

	os.remove(tmp1)
	os.remove(tmp2)
	os.remove(tmp3)
				

inds = get_clusters(c_file, cl)
get_homolog(cl, inds, dir)

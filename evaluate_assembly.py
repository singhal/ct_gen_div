import re
import pandas as pd
import subprocess
import argparse
import os

parser = argparse.ArgumentParser(description="Define homologous loci across species.")
parser.add_argument('--a', help="Assembly for which to run.")
parser.add_argument('--cl', help="Cluster for which to run.")

args = parser.parse_args()
assembly = args.a
cluster = args.cl

dir = '/Volumes/heloderma4/sonal/eco_IBD_oz/snp_calling/%s_alignments/' % (assembly)
outfile = '/Volumes/heloderma4/sonal/eco_IBD_oz/snp_calling/evaluate_assemblies.csv'
# clustering file
c_file = '/Volumes/heloderma4/sonal/eco_IBD_oz/data/clustering/Ctenotus_Lerista.clustering.revised.csv'


def get_clusters(c_file, cluster):
        '''
        get the inds in a given cluster
        '''
        d = pd.read_csv(c_file)
    
        inds = d[d.GMYC_RAxML2 == cluster].sample.tolist()
       
        return inds


def assembly_metrics(dir, cluster, metrics):
	new_dir = dir.replace('alignments', 'assemblies')
	a_file = '%s%s.fa' % (new_dir, cluster)

	metrics['numContigs'] = 0
	metrics['contigLength'] = 0

	f = open(a_file, 'r')
	for l in f:
		if re.search('>', l):
			metrics['numContigs'] += 1
			seq = f.next().rstrip()
			metrics['contigLength'] += len(seq)

	return metrics

	
def mapping_metrics(dir, inds, metrics):
	names = ['MAQ', 'totalReads', 'percentMapped', 'percentPaired', 'nonuniqueReads']
	for name in names:
		metrics[name] = 0

	for ind in inds:
		bam = '%s%s.bwamem.sorted.bam' % (dir, ind)
		p = subprocess.Popen("samtools flagstat %s" % bam, stdout=subprocess.PIPE, shell=True)
		x = [l.rstrip() for l in p.stdout]
		metrics['totalReads'] += int(re.search('^(\d+)', x[0]).group(1))
		metrics['percentMapped'] += float(re.search('([\d\.]+)%', x[4]).group(1))
		metrics['percentPaired'] += float(re.search('([\d\.]+)%', x[8]).group(1))

		p = subprocess.Popen("samtools view %s | awk '$5 < 1' | wc" % bam, stdout=subprocess.PIPE, shell=True)
		x = [l.rstrip() for l in p.stdout]
		metrics['nonuniqueReads'] += int(re.search('\s+(\d+)', x[0]).group(1))

		# mapping quality score excludes non-unique and non mapped
		p = subprocess.Popen("samtools view %s | awk '$5 > 0' | awk '{if ($5 < 255)} {print $5}'" % bam, stdout=subprocess.PIPE, shell=True)
		maq = num_maq = 0
		for l in p.stdout:
			num_maq += 1
			maq += int(l.rstrip())
		metrics['MAQ'] += (maq / float(num_maq))

	for name in names:
		metrics[name] = metrics[name] / float(len(inds)) 
	metrics['percentUnique'] = (metrics['totalReads'] - metrics['nonuniqueReads']) / metrics['totalReads']
	return metrics

def coverage_metrics(dir, inds, metrics):
	bams = ['%s%s.bwamem.sorted.bam' % (dir, ind) for ind in inds]
	out = '%s%s.depth.out' % (dir, cluster)
	subprocess.call("samtools depth %s > %s" % (' '.join(bams), out), shell=True)

	metrics['coveredSites'] = 0
	f = open(out, 'r')
	for l in f:
		d = re.split('\s+', l.rstrip())
		cov = [int(x) for x in d[2:]]
		hi = [x for x in cov if x > 4]
		if (len(hi) / float(len(cov))) >= 0.5:
			metrics['coveredSites'] += 1

	f.close()
	os.remove(out)

	return metrics

metrics = {}
inds = get_clusters(c_file, cluster)
if assembly != 'genome':
	metrics = assembly_metrics(dir, cluster, metrics)
mapping_metrics(dir, inds, metrics)
coverage_metrics(dir, inds, metrics)

o = open(outfile, 'a')
for metric, value in metrics.items():
	o.write('%s,%s,%s,%s\n' % (cluster, assembly, metric, value))
o.close() 

import re
import os
import subprocess
import argparse
import pandas as pd
import gzip

parser = argparse.ArgumentParser(description='Get raw VCF and depth file for clusters.')
parser.add_argument('--cl', help="cluster to run this on")
args = parser.parse_args()
cl = args.cl

c_file = '/Volumes/heloderma4/sonal/eco_IBD_oz/data/clustering/Ctenotus_Lerista.clustering.revised.csv'
seq_dir = '/Volumes/heloderma4/sonal/eco_IBD_oz/snp_calling/rainbow_assemblies/'
bam_dir = '/Volumes/heloderma4/sonal/eco_IBD_oz/snp_calling/rainbow_alignments/'
vcf_dir = '/Volumes/heloderma4/sonal/eco_IBD_oz/snp_calling/variants/'
cov_dir = '/Volumes/heloderma4/sonal/eco_IBD_oz/snp_calling/coverage/'
min_qual = 20
min_cov = 5


def get_clusters(c_file):
        d = pd.read_csv(c_file)
        d = d[d.GMYC_RAxML2.notnull()].groupby('GMYC_RAxML2')

	clusters = dict([(name, sorted(g.sample.tolist()))  for name, g in d])

        return clusters


def get_raw_vcf(vcf_dir, bam_dir, seq_dir, cl, inds):
	raw_vcf = '%s%s.raw.vcf.gz' % (vcf_dir, cl)
	bam_files = ['%s%s.bwamem.final.low_mm.bam' % (bam_dir, ind) for ind in inds]
	seq = '%s%s.fa' % (seq_dir, cl)	

	if not os.path.isfile(raw_vcf):
		subprocess.call("samtools mpileup -ugf %s %s | /Volumes/heloderma4/sonal/bin/bcftools/bcftools call -mO z -o %s" % (seq, ' '.join(bam_files), raw_vcf), shell=True)

	
def get_depth(cov_dir, bam_dir, cl, inds):
	cov_file = '%s%s.depth.out' % (cov_dir, cl)
	bam_files = ['%s%s.bwamem.final.low_mm.bam' % (bam_dir, ind) for ind in inds]

	if not os.path.isfile(cov_file + '.gz'):
		subprocess.call("samtools depth %s > %s" % (' '.join(bam_files), cov_file), shell=True)
		subprocess.call("gzip %s" % cov_file, shell = True)


def initialize_depth(cov_dir):
	out = '%sdepth_profiling.csv' % cov_dir
	if not os.path.isfile(out):
		o = open(out, 'w')
		o.write('cluster,individual,depth,counts\n')
		o.close()
	return out


def summarize_depth(cov_dir, out, cl, inds):
	depthout = '%s%s.depth.out.gz' % (cov_dir, cl)
	f = gzip.open(depthout, 'r')

	depth = {}
	for ind in inds:
		depth[ind] = {}
	
	for l in f:
		d = [int(x) for x in re.split('\s+', l.rstrip())[2:]]
		for ind_d, ind in zip(d, inds):
			if ind_d > 0:
				if ind_d not in depth[ind]:
					depth[ind][ind_d] = 0
				depth[ind][ind_d] += 1
	f.close()

	o = open(out, 'a')
	for ind in depth:
		for d, count in depth[ind].items():
			o.write('%s,%s,%s,%s\n' % (cl, ind, d, count))
	o.close()


def filter_vcf(cl, inds, cov_dir, vcf_dir):
        cov = '%scoverage_summary.csv' % cov_dir
        cov = pd.read_csv(cov)
        cov = dict([(x, y) for x, y in zip(cov.individual, cov.max_coverage)])

	vcfout = '%s%s.raw.vcf.gz' % (vcf_dir, cl)
	vcf = '%s%s.coverage_filtered.vcf.gz' % (vcf_dir, cl)
	depth = '%s%s.depth.out.gz' % (cov_dir, cl)

	ditch = {}
	depth = gzip.open(depth, 'r')
	for l in depth:
		a = re.split('\s+', l.rstrip())
		d = [int(x) for x in a[2:]]
		keep = [True if ind_d >= min_cov and ind_d <= cov[ind] else False for ind_d, ind in zip(d, inds)]
		if sum(keep) != len(keep):
			if a[0] not in ditch:
				ditch[a[0]] = {}
			ditch[a[0]][a[1]] = keep

	depth.close()

	f = gzip.open(vcfout, 'r')
	o = gzip.open(vcf, 'w')

	for l in f:
		if re.search('^#', l):
			o.write(l)
		else:
			d = re.split('\s+', l.rstrip())
			if float(d[5]) >= min_qual:
				keep = [True] * len(d[9:])
				if d[0] in ditch:
					if d[1] in ditch[d[0]]:
						keep = ditch[d[0]][d[1]]
				if sum(keep) > 0:
					for ix, status in enumerate(keep):      
						if not status:
							d[9 + ix] = './.'
					o.write('%s\n' % '\t'.join(d))
                                                         
	f.close()
	o.close()

clusters = get_clusters(c_file)
# get_raw_vcf(vcf_dir, bam_dir, seq_dir, cl, clusters[cl])
# get_depth(cov_dir, bam_dir, cl, clusters[cl])
# out = initialize_depth(cov_dir)
# summarize_depth(cov_dir, out, cl, clusters[cl])
filter_vcf(cl, clusters[cl], cov_dir, vcf_dir)

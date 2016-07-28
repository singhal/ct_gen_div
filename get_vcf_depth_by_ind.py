import re
import os
import subprocess
import argparse
import pandas as pd
import gzip

parser = argparse.ArgumentParser(description='Get raw VCF and depth file for clusters.')
parser.add_argument('--ind', help="ind to run this on")
args = parser.parse_args()
ind = args.ind

c_file = '/Volumes/heloderma4/sonal/eco_IBD_oz/data/clustering/Ctenotus_Lerista.clustering.revised.csv'
seq_dir = '/Volumes/heloderma4/sonal/eco_IBD_oz/snp_calling/rainbow_assemblies/'
bam_dir = '/Volumes/heloderma4/sonal/eco_IBD_oz/snp_calling/rainbow_alignments/'
vcf_dir = '/Volumes/heloderma4/sonal/eco_IBD_oz/snp_calling/ind_variants/'
cov_dir = '/Volumes/heloderma4/sonal/eco_IBD_oz/snp_calling/ind_variants/coverage/'
min_qual = 20
min_cov = 5


def get_clusters(c_file):
        d = pd.read_csv(c_file)
        d = d[d.GMYC_RAxML2.notnull()]

	inds = dict(zip(d.sample, d.GMYC_RAxML2))

        return inds


def get_raw_vcf(vcf_dir, bam_dir, seq_dir, cl, ind):
	raw_vcf = '%s%s.raw.vcf.gz' % (vcf_dir, ind)
	bam_file = '%s%s.bwamem.final.low_mm.bam' % (bam_dir, ind)
	seq = '%s%s.fa' % (seq_dir, cl)	

	if not os.path.isfile(raw_vcf):
		subprocess.call("samtools mpileup -ugf %s %s | /Volumes/heloderma4/sonal/bin/bcftools/bcftools call -mO z -o %s" % (seq, bam_file, raw_vcf), shell=True)

	
def get_depth(cov_dir, bam_dir, cl, ind):
	cov_file = '%s%s.depth.out' % (cov_dir, ind)
	bam_file = '%s%s.bwamem.final.low_mm.bam' % (bam_dir, ind)

	if not os.path.isfile(cov_file + '.gz'):
		subprocess.call("samtools depth %s > %s" % (bam_file, cov_file), shell=True)
		subprocess.call("gzip %s" % cov_file, shell = True)


def initialize_depth(cov_dir):
	out = '%sdepth_profiling.csv' % cov_dir
	if not os.path.isfile(out):
		o = open(out, 'w')
		o.write('cluster,individual,depth,counts\n')
		o.close()
	return out


def summarize_depth(cov_dir, out, cl, ind):
	depthout = '%s%s.depth.out.gz' % (cov_dir, ind)
	f = gzip.open(depthout, 'r')

	depth = {}
	
	for l in f:
		ind_d = int(re.split('\s+', l.rstrip())[2])
		if ind_d > 0:
			if ind_d not in depth:
				depth[ind_d] = 0
			depth[ind_d] += 1
	f.close()

	o = open(out, 'a')
	for d, count in depth.items():
		o.write('%s,%s,%s,%s\n' % (cl, ind, d, count))
	o.close()


def filter_vcf(cl, ind, cov_dir, vcf_dir):
        cov = '%scoverage_summary.csv' % cov_dir
        cov = pd.read_csv(cov)
        cov = dict([(x, y) for x, y in zip(cov.individual, cov.max_coverage)])

	vcfout = '%s%s.raw.vcf.gz' % (vcf_dir, ind)
	vcf = '%s%s.coverage_filtered.vcf.gz' % (vcf_dir, ind)
	depth = '%s%s.depth.out.gz' % (cov_dir, ind)

	ditch = {}
	depth = gzip.open(depth, 'r')
	for l in depth:
		a = re.split('\s+', l.rstrip())
		d = int(a[2])
		keep = False
		if d >= min_cov and d <= cov[ind]:
			keep = True
		if not keep:
			if a[0] not in ditch:
				ditch[a[0]] = {}
			ditch[a[0]][a[1]] = 1

	depth.close()

	f = gzip.open(vcfout, 'r')
	o = gzip.open(vcf, 'w')

	for l in f:
		if re.search('^#', l):
			o.write(l)
		else:
			d = re.split('\s+', l.rstrip())
			if float(d[5]) >= min_qual:
				printIt = True
				if d[0] in ditch:
					if d[1] in ditch[d[0]]:
						printIt = False
				if printIt:		
					o.write('%s\n' % '\t'.join(d))
                                                         
	f.close()
	o.close()

inds = get_clusters(c_file)
get_raw_vcf(vcf_dir, bam_dir, seq_dir, inds[ind], ind)
get_depth(cov_dir, bam_dir, inds[ind], ind)
out = initialize_depth(cov_dir)
# summarize_depth(cov_dir, out, inds[ind], ind)
filter_vcf(inds[ind], ind, cov_dir, vcf_dir)

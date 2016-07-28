import re 
import glob
import pandas as pd
import os
import subprocess
import argparse
import sys

parser = argparse.ArgumentParser(description="Run final alignment steps.")
parser.add_argument('--cl', help="Cluster for which to run script.")
args = parser.parse_args()
cl = args.cl

min_qual = 20

c_file1 = '/Volumes/heloderma4/sonal/eco_IBD_oz/data/clustering/Ctenotus_Lerista.clustering.revised.csv'
d = pd.read_csv(c_file1)
d = d[d.GMYC_RAxML2.notnull()].groupby('GMYC_RAxML2')
clusters = dict([(cluster, sorted(group.sample.tolist())) for cluster, group in d])
inds = clusters[cl]

gatk = '/Volumes/heloderma4/sonal/bin/GenomeAnalysisTK.jar'
dir = '/Volumes/heloderma4/sonal/eco_IBD_oz/snp_calling/rainbow_alignments/'
seq_dir = '/Volumes/heloderma4/sonal/eco_IBD_oz/snp_calling/rainbow_assemblies/'

bamfiles = ['%s%s.bwamem.unique.bam' % (dir, ind) for ind in inds]
for ix, bamfile in enumerate(bamfiles):
	if not os.path.isfile(bamfile):
		new_bamfile = bamfile.replace('unique', 'final')
		if os.path.isfile(new_bamfile):
			bamfiles[ix] = new_bamfile
	
raw_vcf = '%s%s.raw.vcf' % (dir, cl)
filt_vcf = '%s%s.filtered.vcf' % (dir, cl)
seq = '%s%s.fa' % (seq_dir, cl)

subprocess.call("samtools mpileup -ugf %s %s | /Volumes/heloderma4/sonal/bin/bcftools/bcftools call -vmO v -o %s" % (seq, ' '.join(bamfiles), raw_vcf), shell=True)

# filter the vcf
f = open(raw_vcf, 'r')
o = open(filt_vcf, 'w')

for l in f:
	if re.match('^#', l):
		o.write(l)
	else:
		d = re.split('\t', l.rstrip())
		if float(d[5]) >= min_qual:
			o.write(l)

f.close()
o.close()

for orig, ind in zip(bamfiles, inds):
	final = '%s%s.bwamem.final.bam' % (dir, ind)
	if not os.path.isfile(final):
		recal = '%s%s.recal.table' % (dir, ind)

		subprocess.call("samtools index %s" % orig, shell=True)
		subprocess.call('java -Xmx10g -jar %s -T BaseRecalibrator -R %s -knownSites %s -I %s -o %s' % (gatk, seq, filt_vcf, orig, recal), shell=True)
		subprocess.call('java -Xmx10g -jar %s -T PrintReads -R %s -I %s --BQSR %s -o %s' % (gatk, seq, orig, recal, final), shell=True) 

		os.remove(recal)
		os.remove(orig)
os.remove(raw_vcf)
os.remove(filt_vcf)
os.remove(filt_vcf + '.idx')

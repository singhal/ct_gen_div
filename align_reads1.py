import argparse
import os
import subprocess
import pandas as pd

parser = argparse.ArgumentParser(description='reads')
parser.add_argument('--ind', help="ind to run this on")
args = parser.parse_args()
ind = args.ind

c_file = '/Volumes/heloderma4/sonal/eco_IBD_oz/data/clustering/Ctenotus_Lerista.clustering.revised.csv'
seq_dir = '/Volumes/heloderma4/sonal/eco_IBD_oz/snp_calling/rainbow_assemblies/'
read_dir = '/Volumes/heloderma4/sonal/eco_IBD_oz/snp_calling/reads/'
out_dir = '/Volumes/heloderma4/sonal/eco_IBD_oz/snp_calling/rainbow_alignments/'

def get_cluster(c_file, ind):
	d = pd.read_csv(c_file)
	cl = d[d.sample == ind].GMYC_RAxML2.tolist()[0]
	return cl


def prepare_seq(cl, seq_dir):
	seq = '%s%s.fa' % (seq_dir, cl)
	if not os.path.isfile(seq + '.bwt'):
		subprocess.call("/Volumes/heloderma4/sonal/bin/bwa-0.7.12/bwa index %s" % seq, shell=True)
	if not os.path.isfile(seq + '.fai'):
		subprocess.call("samtools faidx %s" % seq, shell=True)
	if not os.path.isfile(seq.replace('.fa', '.dict')):
		subprocess.call("java -jar /Volumes/heloderma4/sonal/bin/picard.jar CreateSequenceDictionary R=%s O=%s" % (seq, seq.replace('.fa', '.dict')), shell=True)


def align_seq(ind, cl, out_dir, read_dir, seq_dir):
	r1 = '%s%s_R1.fastq.gz' % (read_dir, ind)
	r2 = '%s%s_R2.fastq.gz' % (read_dir, ind)
	seq = '%s%s.fa' % (seq_dir, cl)

	out1 = '%s%s.sam' % (out_dir, ind)
	out2 = '%s%s.mateFixed.bam' % (out_dir, ind)
	out3 = '%s%s.mateFixed.sorted.bam' % (out_dir, ind)
	out4 = '%s%s.rg.mateFixed.sorted.bam' % (out_dir, ind)
	intervals = '%s%s.intervals' % (out_dir, ind)
	out5 = '%s%s.realigned.rg.mateFixed.sorted.bam' % (out_dir, ind)
	out6 = '%s%s.bwamem.unique.bam' % (out_dir, ind)

	tmpdir = '%s%s/' % (out_dir, ind)
	if not os.path.isdir(tmpdir):
		os.mkdir(tmpdir)

	# align
	subprocess.call("/Volumes/heloderma4/sonal/bin/bwa-0.7.12/bwa mem -t 4 %s %s %s > %s" % (seq, r1, r2, out1), shell=True)
	# fixmate
	subprocess.call("samtools fixmate -O bam %s %s" % (out1, out2), shell=True)
	# sorted
	subprocess.call("samtools sort -O bam -o %s -T %s %s" % (out3, tmpdir, out2), shell=True)
	# readgroup
	subprocess.call("java -jar /Volumes/heloderma4/sonal/bin/picard.jar AddOrReplaceReadGroups INPUT=%s OUTPUT=%s RGLB=%s RGPL=Illumina RGPU=%s RGSM=%s" % (out3, out4, ind, ind, ind), shell=True)
	subprocess.call("samtools index %s" % out4, shell=True)
	# indeltarget
	subprocess.call("java -Xmx10g -jar /Volumes/heloderma4/sonal/bin/GenomeAnalysisTK.jar -T RealignerTargetCreator -R %s -I %s -o %s -nt 4" % (seq, out4, intervals), shell=True)
	# indelrealigner
	subprocess.call("java -Xmx10g -jar /Volumes/heloderma4/sonal/bin/GenomeAnalysisTK.jar -T IndelRealigner -R %s -I %s -targetIntervals %s -o %s" % (seq, out4, intervals, out5), shell=True)
	# unique
	subprocess.call("samtools view -b -q 1 %s > %s" % (out5, out6), shell=True)

	# call = [os.remove(x) for x in [out1, out2, out3, out4, out5, intervals, out4 + '.bai', out5.replace('bam', 'bai')]]
	os.rmdir(tmpdir) 

# get cluster
cl = get_cluster(c_file, ind)
# prepare seqfiles
prepare_seq(cl, seq_dir)
# align it all the way until time to call SNPs
align_seq(ind, cl, out_dir, read_dir, seq_dir)

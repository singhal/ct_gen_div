import glob
import subprocess
import argparse

parser = argparse.ArgumentParser(description='reads')
parser.add_argument('--ind', help="ind to run this on")
args = parser.parse_args()
ind = args.ind

dir = '/Volumes/heloderma4/sonal/eco_IBD_oz/snp_calling/reads/'
outdir = '/Volumes/heloderma4/sonal/eco_IBD_oz/snp_calling/trimmed_reads/'

def trim_reads(dir, outdir, file1):
	file2 = file1.replace('_R1', '_R2')
	out1 = file1.replace(dir, outdir)
	out2 = file2.replace(dir, outdir)
	outu1 = out1.replace('_R1', '_u1')
	outu2 = out2.replace('_R2', '_u2')

	subprocess.call('java -jar /Volumes/heloderma4/sonal/bin/Trimmomatic-0.33/trimmomatic-0.33.jar PE -threads 12 -phred33 %s %s %s %s %s %s SLIDINGWINDOW:4:25 MINLEN:36' % (file1, file2, out1, outu1, out2, outu2), shell=True)

file = '%s%s_R1.fastq.gz' % (dir, ind)
trim_reads(dir, outdir, file)

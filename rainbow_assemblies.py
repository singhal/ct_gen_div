import re
import glob
import argparse
import os
import subprocess

parser = argparse.ArgumentParser(description='reads')
parser.add_argument('--ind', help="ind to run this on")
args = parser.parse_args()
ind = args.ind

out_dir = '/Volumes/heloderma4/sonal/eco_IBD_oz/snp_calling/rainbow_assemblies/indiv/'

file1 = '/Volumes/heloderma4/sonal/eco_IBD_oz/snp_calling/trimmed_reads/%s_R1.fastq.gz' % ind
file2 = file1.replace('R1', 'R2')
	
out1 = '%s%s.rainbow1' % (out_dir, ind)
out2 = '%s%s.rainbow2' % (out_dir, ind)
out3 = '%s%s.rainbow3' % (out_dir, ind)
out4 = '%s%s.fa' % (out_dir, ind)	
	
subprocess.call('/Volumes/heloderma4/sonal/bin/rainbow_2.0.4/rainbow cluster -1 %s -2 %s > %s\n' % (file1, file2, out1), shell=True)
subprocess.call('/Volumes/heloderma4/sonal/bin/rainbow_2.0.4/rainbow div -i %s -o %s\n' % (out1, out2), shell=True)
subprocess.call('/Volumes/heloderma4/sonal/bin/rainbow_2.0.4/rainbow merge -i %s -a -o %s\n' % (out2, out3), shell=True)
subprocess.call('/Volumes/heloderma4/sonal/bin/rainbow_2.0.4/select_best_rbcontig_plus_read1.pl %s %s > %s\n' % (out3, out2, out4), shell=True)

os.remove(out1)
os.remove(out2)
os.remove(out3)


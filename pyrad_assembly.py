import re
import gzip
import subprocess
import argparse
import os
import itertools
import shutil

parser = argparse.ArgumentParser(description="Run PyRAD assembly.")
parser.add_argument('--ind', help="Individual for which to run assembly.")
args = parser.parse_args()
ind = args.ind

outdir = '/Volumes/heloderma4/sonal/eco_IBD_oz/snp_calling/pyrad_assemblies/'
readdir = '/Volumes/heloderma4/sonal/eco_IBD_oz/snp_calling/trimmed_reads/'

def prepare_pyrad(outdir, ind):
	inddir = '%s%s/'% (outdir, ind)
	if not os.path.exists(inddir):
		os.makedirs(inddir)
	editdir = '%sedits/' % inddir
	if not os.path.exists(editdir):
		os.makedirs(editdir)

	os.chdir(inddir)
	subprocess.call("pyRAD -n", shell=True)

	subprocess.call("sed -i '' -e 's|.*## 4.|/Volumes/heloderma4/sonal/bin/vsearch/vsearch   ## 4.|' params.txt", shell=True)
	subprocess.call("sed -i '' -e 's/.*## 10./0.98   ##10./' params.txt", shell=True)
	subprocess.call("sed -i '' -e 's/.*## 11./pairddrad   ## 11./' params.txt", shell=True)
	subprocess.call("sed -i '' -e 's|.*## 18.|%s   ## 18.|' params.txt" % inddir, shell=True)

	return (inddir, editdir)

def run_pyrad(inddir):
	os.chdir(inddir)
	subprocess.call("pyRAD -p params.txt -s 3", shell=True)

def make_starting_file(ind, readdir, inddir, editdir):
	read1 = '%s%s_R1.fastq.gz' % (readdir, ind)
	read2 = '%s%s_R2.fastq.gz' % (readdir, ind)
	out = '%s%s.edit' % (editdir, ind)

	f1 = gzip.open(read1, 'r')
	f2 = gzip.open(read2, 'r')
	o = open(out, 'w')

	ix = 0
	for l1 in f1:
		l1 = l1.rstrip()
		l2 = f2.next().rstrip()

		if re.search('^@', l1):
			o.write('>%s_pair_%s\n' % (ind, ix))
	
			l1 = f1.next().rstrip()
			l2 = f2.next().rstrip() 
			o.write('%snnnn%s\n' % (l1, l2))

			l1 = f1.next().rstrip()
        	        l2 = f2.next().rstrip() 

			l1 = f1.next().rstrip()
        	        l2 = f2.next().rstrip()
			ix += 1

	f1.close()
	f2.close()
	o.close()

def make_final_file(inddir, outdir, ind):
	file_in = '%sclust0.98/%s.clustS.gz' % (inddir, ind)
	file_out = '%s%s.fa' % (outdir, ind)

	f = gzip.open(file_in, 'r')
	o = open(file_out, 'w')

	tmp_loci = {}
	for l in f:
		if re.search('>', l):
			seq = f.next().rstrip()
			vals = re.split(';', re.search('>(\S+)', l).group(1))
			tmp_loci[vals[0]] = {'seq': seq, 'size': vals[1]}
		elif re.search('//', l):
			if tmp_loci:
				loci = sorted(tmp_loci, key=lambda x: tmp_loci[x]['size'], reverse=True)
				seq = re.sub('-', '', tmp_loci[loci[0]]['seq'])
				o.write('>%s\n%s\n' % (loci[0], seq))
			tmp_loci = {}
		else: 
			print "ERROR! Something is wrong with clustS file."
	f.close()
	o.close()

	shutil.rmtree(inddir)

inddir, editdir = prepare_pyrad(outdir, ind)
make_starting_file(ind, readdir, inddir, editdir)
run_pyrad(inddir)
make_final_file(inddir, outdir, ind)

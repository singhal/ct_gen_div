import os
import re
import pandas as pd
import subprocess
import glob
import argparse

parser = argparse.ArgumentParser(description="Make gene tree files.")
parser.add_argument('--g', help="Genus for which to run.")
parser.add_argument('--m', help="Missing level tolerated.")
args = parser.parse_args()

genus = args.g
missing = float(args.m)

outdir = '/Volumes/heloderma4/sonal/eco_IBD_oz/data/species_tree/%s/' % (genus)
seqdir = '/Volumes/heloderma4/sonal/eco_IBD_oz/snp_calling/rainbow_assemblies/'

# homology file
hom_file = '%s%s_homology_across_species.txt' % (seqdir, genus)

# make dir if it doesn't exist
if not os.path.exists(outdir):
	os.makedirs(outdir)


def define_loci(hom_file, seqdir, genus, missing):
	files = glob.glob('%s%s*fa' % (seqdir, genus))
	inds = [re.search('([a-z|A-Z|_|0-9]+).fa', file).group(1) for file in files]

	num_ind = len(inds)

	d = pd.read_csv(hom_file, sep = '\t')	
	d = d[d.numMatches >= (num_ind * missing)]
	
	loci = {}
	for c, match in zip(d.contig, d.matches):
		loci[c] = {}
		matches = re.split(',', match)
		matches = [re.split(':', match) for match in matches]
		for match in matches:
			loci[c][match[0]] = match[1]

	return loci, files, inds


def get_sequences(files):
	seqs = {}

	for file in files:
		f = open(file, 'r')
		for l in f:
			if re.search('>', l):
				id = re.search('>(\S+)', l).group(1)
				seq = f.next().rstrip()
				seqs[id] = seq
		f.close()

	return seqs	


def revcomp(seq):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
	rc = "".join(complement.get(base) for base in reversed(seq))
	return rc


def make_alignments(outdir, loci, seqs, inds):
	summary_file = '%sloci_summary.csv' % outdir
	s = open(summary_file, 'w')
	s.write('locus,sequenced\n')

	for locus in loci:
		out = '%s%s.fa' % (outdir, locus)
		alnout = '%s%s.fa.aln' % (outdir, locus)

		s.write('%s,%.3f\n' % (locus, len(loci[locus]) / float(len(inds))))

		o = open(out, 'w')
		for contig, strand in loci[locus].items():
			seq = seqs[contig]
			if strand == '-':
				seq = revcomp(seq)
			species = re.search('^(\S+)_\d+$', contig).group(1)
			o.write('>%s\n%s\n' % (species, seq))
		o.close()

		subprocess.call("muscle -in %s -out %s -quiet" % (out, alnout), shell=True)
		# subprocess.call("mafft --adjustdirection --auto --quiet %s > %s" % (out, alnout), shell=True)
		os.remove(out)
	s.close()


loci, files, inds = define_loci(hom_file, seqdir, genus, missing)
seqs = get_sequences(files)
make_alignments(outdir, loci, seqs, inds)

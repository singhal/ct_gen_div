import os
import re
import pandas as pd
import subprocess
import glob
import argparse
import random

parser = argparse.ArgumentParser(description="Run gene tree analysis.")
parser.add_argument('--l', help="Locus for which to run.")
parser.add_argument('--g', help="Genus for which to run.")
args = parser.parse_args()

genus = args.g
locus = args.l

def convert_phyml(locus_file):
	f = open(locus_file, 'r')
	phy_file = locus_file.replace('.fa.aln', '.aln.phy')
	o = open(phy_file, 'w')

	seq = {}
	id = ''
	for l in f:
		if re.search('>', l):
			id = re.search('>(\S+)', l.rstrip()).group(1)
			seq[id] = ''
		else:
			seq[id] += l.rstrip()
	f.close()

	o.write(' %s %s\n' % (len(seq), len(seq.values()[0])))
	for sp, s in seq.items():
		o.write('%s   %s\n' % (sp, s))
	o.close()

	return phy_file


def run_raxml(locus, outdir, phy_file):
	os.chdir(outdir)
	subprocess.call('/Volumes/heloderma4/sonal/bin/RAxML/raxmlHPC -x %s -# 100 -p %s -m GTRCAT -f a -n %s -s %s' % (random.randint(0,1000),  random.randint(0,1000), locus, phy_file), shell=True)

	orig_boot = 'RAxML_bootstrap.%s' % locus
	orig_tree = 'RAxML_bipartitions.%s' % locus

	new_boot = '%s.bootstrap.trees' % locus
	new_tree = '%s.bestTree.tre' % locus

	os.rename(orig_boot, new_boot)
	os.rename(orig_tree, new_tree)

	subprocess.call("rm RAxML_*%s" % locus, shell=True) 


dir = '/Volumes/heloderma4/sonal/eco_IBD_oz/data/species_tree/%s/' % genus
locus_file = '%s%s.fa.aln' % (dir, locus)
outdir = '/Volumes/heloderma4/sonal/eco_IBD_oz/data/species_tree/gene_trees_%s/' % genus

phy_file = convert_phyml(locus_file)
run_raxml(locus, outdir, phy_file)

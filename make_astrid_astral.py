import re
import pandas as pd
import os

# genus = ['Le', 'Ct']
genera = ['Ct']
missing = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3]
ALIGN_QUAL = 0.7
astral = '/Volumes/heloderma4/sonal/bin/ASTRAL/Astral/astral.4.7.8.jar'

dir = '/Volumes/heloderma4/sonal/eco_IBD_oz/data/species_tree/'
outdir = '%standy/' % dir

def get_loci(genus, miss):
	file = '%s/%s/loci_summary.csv' % (dir, genus)
	d = pd.read_csv(file)
	loci = d[(d.sequenced >= miss) & (d.alignment >= ALIGN_QUAL)].locus.tolist()
	return loci


def make_conc_trees(genus, miss, loci, outdir):
	outfile = '%s%s_missing%s.bestML_poly.trees' % (outdir, genus, miss)
	o = open(outfile, 'w')
	for locus in loci:
		tree = '%sgene_trees_%s/%s.bestTree_poly.tre'  % (dir, genus, locus)
		if os.path.isfile(tree):
			f = open(tree, 'r')
			o.write(f.next())
			f.close()
	o.close()
	return outfile


def make_astral_bs(genus, miss, loci, outdir):
	outfile = '%sastral/%s_missing%s.bs_poly.txt' % (outdir, genus, miss)
	o = open(outfile, 'w')
        for locus in loci:
                tree = '%sgene_trees_%s/%s.bootstrap.trees' % (dir, genus, locus)
                if os.path.isfile(tree):
			tree = '/home/sosi/astral/gene_trees_%s/%s.bootstrap.trees' % (genus, locus)
                        o.write(tree + '\n')
        o.close()
        return outfile


def make_astrid_bs(genus, miss, loci, outdir):
        outfiles = []
	bs = {}
        for locus in loci:
                tree = '%sgene_trees_%s/%s.bootstrap.trees' % (dir, genus, locus)
                if os.path.isfile(tree):
			t = open(tree, 'r')
			bs[locus] = [x.rstrip() for x in t.readlines()]
			t.close()
	for i in range(0, 100):
		outfile = '%sastrid/%s_missing%s.bs_poly%s' % (outdir, genus, miss, i)
	        outfiles.append(outfile)
		o = open(outfile, 'w')
		for locus in loci:
			if locus in bs:
				o.write(bs[locus][i] + '\n')
        	o.close()
        return outfiles


def astral_out(treefile, astral_bsfile, genus, missing):
	outfile = '%sastral.sh' % outdir
	o = open(outfile, 'a')
	out = '%sastral/%s_missing%s' % (outdir, genus, missing)
	o.write("java -Xmx4000M -jar %s -i %s -b %s -r 100 -o %s -k completed\n" % (astral, treefile, astral_bsfile, out))
	o.close()


def astrid_out(treefile, astrid_bsfiles, genus, missing):
        outfile = '%sastrid.sh' % outdir
        o = open(outfile, 'a')
        out = '%sastrid/%s_missing%s' % (outdir, genus, missing)
        o.write("ASTRID -i %s -o %s -m bionj\n" % (treefile, out))
	for i, bsfile in enumerate(astrid_bsfiles):
		out = '%sastrid/%s_missing%s_%s' % (outdir, genus, missing, i)
        	o.write("ASTRID -i %s -o %s -m bionj\n" % (bsfile, out))
        o.close()


for genus in genera:
	for miss in missing:
		loci = get_loci(genus, miss)
		# treefile = make_conc_trees(genus, miss, loci, outdir)
		astral_bsfile = make_astral_bs(genus, miss, loci, outdir)
		# astrid_bsfiles = make_astrid_bs(genus, miss, loci, outdir)
		# astral_out(treefile, astral_bsfile, genus, miss)
		# astrid_out(treefile, astrid_bsfiles, genus, miss)

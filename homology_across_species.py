import re
import subprocess
import os
import argparse
import random
import string
import sys
import pandas as pd

parser = argparse.ArgumentParser(description="Define homologous loci across species.")
parser.add_argument('--dir', help="Directory for which to run.")
parser.add_argument('--genus', help="Genus for which to run.")

args = parser.parse_args()
dir = args.dir
genus = args.genus

c_file = '/Volumes/heloderma4/sonal/eco_IBD_oz/data/clustering/Ctenotus_Lerista.clustering.revised.csv'
d = pd.read_csv(c_file)
d = d[d.GMYC_RAxML2.notnull()]
clusters = d[d.GMYC_RAxML2.str.contains(genus)].GMYC_RAxML2.unique()

files = ['%s%s.fa' % (dir, cl) for cl in clusters]
WCLUST = 0.8

def create_starter(dir, file, genus, ix):
	homhash = {}

	starting = '%s%s.tmp.fa' % (dir, genus)
	f = open(file, 'r')
	o = open(starting, 'w')
	for l in f:
		if re.search('>', l):
			id = re.search('>(\S+)', l).group(1)
			newid = '%s_%s' % (genus, ix)

			homhash[newid] = {}
			homhash[newid][id] = '+'
			
			seq = f.next().rstrip()
			o.write('>%s\n%s\n' % (newid, seq))
			ix += 1
	f.close()
	o.close()

	return starting, homhash, ix


def vsearch(dir, tmpfile, file, genus, num):
	out = '%s%s_%s_search' % (dir, genus, num)
	subprocess.call("vsearch --usearch_global %s --db %s --userout %s --id %s --userfields query+target+evalue+id+qstrand --strand both --threads 4" % (file, tmpfile, out, WCLUST), shell=True)

	return out


def create_new_tmp(dir, tmpfile, file, results, homhash, genus, ix):
	matches1 = {}
	matches2 = {}

	f = open(results, 'r')
	for l in f:
		d = re.split('\s+', l.rstrip())
		# is this step necessary?
		# makes sure it is 1 to 1
		if d[1] not in matches1 and d[0] not in matches2:
			matches1[d[1]] = {'match': d[0], 'perc': float(d[3]), 'strand': d[4]}
			matches2[d[0]] = {'match': d[1], 'perc': float(d[3]), 'strand': d[4]}
		elif d[0] in matches2 and d[1] not in matches1:
			if float(d[3]) > matches2[d[0]]['perc']:
                        	matches1[d[1]] = {'match': d[0], 'perc': float(d[3]), 'strand': d[4]}
                        	matches2[d[0]] = {'match': d[1], 'perc': float(d[3]), 'strand': d[4]}
	f.close()
	os.remove(results)

	for c in matches2:
		homhash[matches2[c]['match']][c] = matches2[c]['strand']

	f = open(file, 'r')
	o = open(tmpfile, 'a')
	for l in f:
		if re.search('>', l):
			id = re.search('>(\S+)', l.rstrip()).group(1)
			seq = f.next().rstrip()
			if id not in matches2:
				new_id = '%s_%s' % (genus, ix)
				ix += 1

				homhash[new_id] = {}	
				homhash[new_id][id] = '+'
		
				o.write('>%s\n%s\n' % (new_id, seq))
	f.close()
	o.close()

	return (tmpfile, homhash, ix)

ix = 0
tmpfile, homhash, ix = create_starter(dir, files[0], genus, ix)
for num, file in enumerate(files[1:]):
	results = vsearch(dir, tmpfile, file, genus, num)
	(tmpfile, homhash, ix) = create_new_tmp(dir, tmpfile, file, results, homhash, genus, ix)
os.remove(tmpfile)

o = open('%s%s_homology_across_species.txt' % (dir, genus), 'w')
o.write('contig\tmatches\tnumMatches\n')
for c, matches in homhash.items():
	matches = ['%s:%s' % (match, homhash[c][match]) for match in matches]
	o.write('%s\t%s\t%s\n' % (c, ','.join(matches), len(matches)))
o.close()



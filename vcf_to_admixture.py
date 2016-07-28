import argparse
import re
import pandas as pd
import gzip

parser = argparse.ArgumentParser(description='Cluster for which to create ADMIXTURE file.')
parser.add_argument('--cl', help="cluster to run this on")
args = parser.parse_args()
cl = args.cl

vcf_dir = '/Volumes/heloderma4/sonal/eco_IBD_oz/snp_calling/variants/'
out_dir = '/Volumes/heloderma4/sonal/eco_IBD_oz/data/admixture/'

recode = {'0/0': 0, '0/1': 1, '1/1': 2, './.': 9}

# get the vcf file
vcf = '%s%s.final.vcf.gz' % (vcf_dir, cl)
f = gzip.open(vcf, 'r')
snps = {}
for l in f:
	if re.search('^#CHROM', l):
		inds = re.split('\t', l.rstrip())[9:]
	elif not re.search('^#', l):
		d = re.split('\t', l.rstrip())
		if d[4] in ['A', 'T', 'C', 'G']:
			genos = [re.search('^(\S\S\S)', x).group(1) for x in d[9:]]
			genos = [recode[x] for x in genos]

			if (genos.count(0) + genos.count(1)) > 0 and (genos.count(1) + genos.count(2)) > 0:
				# only keep the snp if missing data is less than a third
				if (genos.count(9) / float(len(inds))) < 0.33:
					if d[0] not in snps:
						snps[d[0]] = {}
					snps[d[0]][d[1]] = genos
f.close()

winners = {}
for contig in snps:
	loci = {}
	for snp in snps[contig]:
		# count the number of missing sites for a given snp
		loci[snp] = snps[contig][snp].count(9)
	winners[contig] = min(loci, key=loci.get)
contigs = sorted(winners.keys())

outfile = '%s%s.geno' % (out_dir, cl)
o = open(outfile, 'w')
for c in snps:
	geno = ''.join([str(x) for x in snps[c][winners[c]]])
	o.write(geno + '\n')
o.close()

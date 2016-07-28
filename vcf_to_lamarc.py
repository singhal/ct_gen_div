import argparse
import gzip
import re
import random

parser = argparse.ArgumentParser(description='cl')
parser.add_argument('--cl', help="cluster to run this on")
args = parser.parse_args()
cl = args.cl

# min sites; if individual isn't represented at this many sites, drop
min_sites = 1e6
# for a site
per_missing = 0.3

# key files & dirs
vcfdir = '/Volumes/heloderma4/sonal/eco_IBD_oz/snp_calling/variants/'
outdir = '/Volumes/heloderma4/sonal/eco_IBD_oz/data/lamarc/'
base_xml = '%sbase.xml' % outdir

def calc_missing(vcf, min_sites):
	f = gzip.open(vcf)
	for l in f:
		if re.search('#CHROM', l):
			inds = re.split('\t', l.rstrip())[9:]
			inds = dict([(ind, 0) for ind in inds])
		elif not re.search('#', l):
			d = re.split('\t', l.rstrip())
			for geno, ind in zip(d[9:], inds):
				if not re.search('\./\.', geno):
					inds[ind] += 1
	f.close()

	to_drop = [ind for ind in inds if inds[ind] < min_sites]
	inds = sorted(inds.keys())
	
	keep = [ix + 9 for ix, ind in enumerate(inds) if ind not in to_drop]
	inds = [ind for ind in inds if ind not in to_drop]	

	return inds, keep


def get_snps(vcf, keep):
	snps = []
	length = 0

	f = gzip.open(vcf)
	for l in f:
		if not re.search('#', l):
			d = re.split('\t', l.rstrip())
			genos = [d[ix] for ix in keep]
		
			# is missing too high?
			miss = sum([1 for geno in genos if re.search('\./\.', geno)])
			miss = miss / float(len(genos))

			if miss < per_missing:
				length += 1

			if miss < per_missing and d[4] != '.':
				alleles = [d[3]] +  re.split(',', d[4])
				alleles = dict([(str(ix), allele) for ix, allele in enumerate(alleles)])				
				alleles['.'] = 'N'

				# remove indels
				indel = False
				for allele in alleles.values():
					if len(allele) > 1:
						indel = True

				if not indel:
					s = []
					for geno in genos:
						geno = re.search('(\S/\S)', geno).group(1)
						geno = re.split('/', geno)
						random.shuffle(geno)
						geno = [alleles[g] for g in geno]
						s += geno
					snps.append(s)
	f.close()
	return snps, length


def print_snps(inds, snps):
	# turn into haplo
	inds2 = []
	for ix, ind in enumerate(inds):
		inds2.append('a%s_%s' % (ix, ind))
		inds2.append('b%s_%s' % (ix, ind))
	inds2 = [ind[0:10] for ind in inds2]

	seq = {}
	for ix, ind in enumerate(inds2):
		seq[ind] = ''.join([pos[ix] for pos in snps])

	return seq


def print_xml(seq, length, cl, outdir, base):
	out = '%s%s_infile.xml' % (outdir, cl)

	base = open(base, 'r')
	base = base.readlines()

	for ix, line in enumerate(base):
		if re.search("XXXX", line):
			base[ix] = re.sub('XXXX', cl, line)
		elif re.search('LENGTH', line):
			base[ix] = re.sub('LENGTH', str(length), line)

	o = open(out, 'w')
	call = [o.write(x) for x in base]

	for ind, s in seq.items():
		o.write('      <individual name="%s">\n' % ind)
		o.write('       <sample name="%s_0">\n' % ind)
		o.write('        <datablock type="SNP">\n')
		o.write('         %s\n' % s)
		o.write('        </datablock>\n')
		o.write('       </sample>\n')
		o.write('      </individual>\n')
	o.write('     </population>\n')
	o.write('    </region>\n')
	o.write('   </data>\n')
	o.write('  </lamarc>\n')
	o.close()


vcf = '%s%s.final.vcf.gz' % (vcfdir, cl)
inds, keep = calc_missing(vcf, min_sites)
snps, length = get_snps(vcf, keep)
seq = print_snps(inds, snps)
print_xml(seq, length, cl, outdir, base_xml)

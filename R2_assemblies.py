import re
import glob
import pandas as pd
import subprocess
import os

# clustering file
c_file = '/Volumes/heloderma4/sonal/eco_IBD_oz/data/clustering/Ctenotus_Lerista.clustering.revised.csv'
# file with trimmed reads
read_dir = '/Volumes/heloderma4/sonal/eco_IBD_oz/snp_calling/trimmed_reads/'
out_dir = '/Volumes/heloderma4/sonal/eco_IBD_oz/snp_calling/R2_assemblies/'

# how deep do we need the locus to be?
DEPTH = 3
WCLUST = 0.95

golden = [ 'SAMAR_29946_Ct_dux', 'NA_ABTC12574_Le_dese', 'SAMAR_42788_Ct_hebe', 'QMJ_62380_Ct_dux', 
                        'CUMV_14637_Le_dese', 'SAMAR_42343_Ct_stra', 'SAMR_65347_Ct_hebe', 'WAMR_166720_Le_dese', 
                        'SAMR_45029_Ct_stra', 'CUMV_14673_Le_dese', 'CUMV_14635_Le_dese', 'CUMV_14362_Ct_dux', 
                        'CUMV_14667_Le_dese', 'SAMAR_44726_Ct_sept', 'CUMV_14671_Le_dese', 'SAMAR_34780_Ct_stra', 
                        'WAMR_140416_Le_dese', 'CUMV_14369_Ct_dux', 'WAMR_135139_Le_dese', 'SAMAR_42843_Ct_asta', 
                        'CUMV_14616_Le_dese', 'SAMR_46165_Ct_dux', 'CUMV_14618_Le_dese', 'SAMR_55734_Ct_stra', 
                        'CUMV_14619_Le_dese', 'NA_ABTC113897_Ct_hebe', 'QM_84335_Ct_quin', 'SAMR_44903_Ct_stra', 
                        'SAMAR_50902_Ct_stra', 'WAMR_117169_Le_punc', 'SAMR_65380_Ct_hebe', 'CUMV_14672_Le_dese', 
                        'SAMR_36098_Ct_dux', 'WAMR_156155_Le_macr', 'SAMAR_42823_Ct_asta', 'WAMR_117239_Le_punc', 
                        'SAMAR_37916_Ct_stra', 'WAMR_145388_Le_macr', 'SAMAR_46517_Le_dese', 'CUMV_14357_Ct_dux', 
                        'SAMR_42801_Ct_asta', 'SAMAR_42955_Ct_stra', 'SAMR_55678_Ct_hebe', 'SAMAR_42840_Ct_asta', 
                        'QMJ_62373_Ct_asta', 'WAMR_145938_Le_punc', 'SAMAR_42903_Ct_asta', 'WAMR_172265_Le_dese', 
                        'NTMR_22166_Ct_robu', 'SAMR_42487_Ct_stra', 'SAMAR_35906_Le_dese', 'CUMV_14617_Le_dese', 
                        'WAMR_166441_Ct_dux', 'CUMV_14374_Ct_dux', 'SAMAR_46254_Le_dese', 'CUMV_14613_Le_dese', 
                        'CUMV_14620_Le_dese', 'SAMR_45608_Ct_dux', 'SAMAR_42743_Ct_stra', 'SAMR_44387_Ct_dux', 
                        'SAMAR_48573_Ct_dux', 'CUMV_14666_Le_dese', 'CUMV_14365_Ct_dux', 'CUMV_14664_Le_dese', 
                        'SAMAR_42776_Ct_hebe', 'WAMR_163385_Le_dese', 'SAMR_44955_Ct_stra', 'SAMR_56476_Ct_dux', 
                        'SAMAR_42800_Ct_stra', 'CUMV_14366_Ct_dux', 'SAMR_65346_Ct_hebe', 'CUMV_14665_Le_dese', 
                        'CUMV_14638_Le_dese', 'WAMR_116872_Le_hump', 'CUMV_14614_Le_dese', 'WAMR_145352_Le_dese', 
                        'SAMR_62031_Ct_dux', 'CUMV_14669_Le_dese', 'CUMV_14373_Ct_dux', 'CUMV_14674_Le_dese', 
                        'WAMR_102274_Le_macr', 'SAMAR_55677_Ct_hebe', 'SAMR_46857_Ct_stra', 'SAMR_46967_Ct_sept', 
                        'CUMV_14370_Ct_dux', 'CUMV_14639_Le_dese', 'SAMAR_32100_Le_dese', 'CUMV_14670_Le_dese', 
                        'SAMR_65379_Ct_hebe', 'SAMAR_42880_Ct_hebe', 'WAMR_117238_Le_dese', 'WAMR_120847_Le_hump', 
                        'SAMR_42034_Ct_dux', 'CUMV_14668_Le_dese', 'CUMV_14636_Le_dese', 'CUMV_14364_Ct_dux' ]


def get_clusters(c_file, golden):
        '''
        get the inds in a given cluster
        '''
        d = pd.read_csv(c_file)
        
	d = d[d.sample.isin(golden)] 

	# group by relevant cluster
        group = d.groupby('GMYC_RAxML2')
       
        # some clusters will only have
        # one individual in them
        clusters = dict([(name, d.sample.tolist()) for name, d in group])
        return clusters


def make_read_files(clusters, read_dir, out_dir):
	for cluster, inds in clusters.items():
		out1 = '%s%s.tmp1.fasta' % (out_dir, cluster)
		out2 = '%s%s.fasta' % (out_dir, cluster)
		if not (os.path.isfile(out1) or os.path.isfile(out2)):
			for ind in inds:
				file1 = '%s%s_R1.fastq.gz' % (read_dir, ind)
				file2 = file1.replace('_R1', '_R2')

				subprocess.call("pear -f %s -r %s -o %s%s -j 4\n" % (file1, file2, out_dir, ind), shell=True)
			finalq = '%s%s.tmp1.fastq' % (out_dir,cluster)
			final = '%s%s.tmp1.fasta' % (out_dir, cluster)	
			subprocess.call("cat %s*.assembled.fastq %s*.unassembled.reverse.fastq > %s" % (out_dir, out_dir, finalq), shell=True)
			subprocess.call("seqtk seq -A -L 60 %s > %s" % (finalq, final), shell=True)
			subprocess.call("rm %s*fastq" % out_dir, shell=True)


def derep_cluster_file(clusters, out_dir):
	for cluster, inds in clusters.items():
		start = '%s%s.tmp1.fasta' % (out_dir, cluster)
		middle1 = '%s%s.tmp2.fasta' % (out_dir, cluster)
		middle2 = '%s%s.tmp3.fasta' % (out_dir, cluster)
		final = '%s%s.fasta' % (out_dir, cluster)

		if not os.path.isfile(final):
			subprocess.call("/Volumes/heloderma4/sonal/bin/vsearch/vsearch --derep_fulllength %s --output %s --sizeout --threads 4" % (start, middle1), shell=True)
			subprocess.call("/Volumes/heloderma4/sonal/bin/vsearch/vsearch --cluster_fast %s --centroids %s --id %s --sizeout --sizein --threads 4" % (middle1, middle2, WCLUST), shell=True)

			cluster_depth = len(inds) * DEPTH
			if cluster_depth > 10:
				cluster_depth = 10

			f = open(middle2, 'r')
			o = open(final, 'w')
			for l in f:
				if re.search('>', l):
					seq = f.next().rstrip()
					size = re.search('size=(\d+)', l).group(1)
	                		if int(size) >= cluster_depth:
	                        		o.write(l.rstrip() + '\n')
						o.write(seq + '\n')
			f.close()
			o.close()
	
			os.remove(start)
			os.remove(middle1)
			os.remove(middle2)


clusters = get_clusters(c_file, golden)
make_read_files(clusters, read_dir, out_dir)
derep_cluster_file(clusters, out_dir)

import sys,os
import urllib.request 
url_prefix = 'https://pmn.plantcyc.org/gene?orgid=ARA&id='

Gene = open('/mnt/home/peipeiw/Documents/Homo_Chr_Gene_cluster/Ath/Expression/Ath_transcript_length_Araport11.txt','r').readlines()
D = {} 
for inl in Gene[1:]:
	gene = inl.split('\t')[0]
	url = url_prefix + gene
	try:
		urllib.request.urlretrieve(url, "/mnt/home/peipeiw/Documents/Homo_Chr_Gene_cluster/Ath/Pathway_annotation/Try_using_url_search_PMN/%s_pathway_info.txt"%gene)
		if gene not in D:
			D[gene] = {}
			D[gene]['pathway'] = []
			D[gene]['evidence'] = []
		inp = open("/mnt/home/peipeiw/Documents/Homo_Chr_Gene_cluster/Ath/Pathway_annotation/Try_using_url_search_PMN/%s_pathway_info.txt"%gene,'r').readlines()
		for inl in inp:
			if "PATHWAY&object=" in inl:
				pathway = inl.split('PATHWAY&object=')[1].split('\"')[0]
				D[gene]['pathway'].append(pathway)
			if "EV-EXP" in inl:
				tem = inl.split('EV-EXP')
				for t in tem[1:]:
					if '.gif' not in t:
						evi = 'EV-EXP' + t.split(':')[0]
						D[gene]['evidence'].append(evi)
	except:
		print('Errors for %s\n'%gene)

out = open('Ath_gene_pathway_evidence.txt','w')
for gene in D:
	evidences = ';'.join(D[gene]['evidence'])
	for pathway in D[gene]['pathway']:
		out.write('%s\t%s\t%s\n'%(gene,pathway,evidences))

out.close()

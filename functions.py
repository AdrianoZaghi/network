from classes import datab
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt


def grow_tax_level(in_datab, out_level, taxonomy_file):
	"""This function returns a new datab object where the occurrence is reorganized in therms of an higher level classification
		in_datab 	 	 	 	 	 the original datab object from which will be taken the occurrence
		out_level 	 	 	 	 	 the higher level of classification required
		taxonomy_file 	 	 	 	 a csv file where is specified the relations between different level classes
	The only usefull values for out_level are those in the column labels of the taxonomy_file
	This file must contain also a label matching with the "Classification level" specified in in_datab.description
	Moreove the classification levels in the taxonomy_file must be ordered from right to the left for increasing granularity
	The out level can't be at the righ of the classification level specified in in_datab.description
	"""
	with open(taxonomy_file, mode = 'r', encoding = "utf-8") as r:
		ref = pd.read_csv(r, index_col=0)
		if out_level not in ref.columns:
			print("Classification required not contemplated in the file of reference")
			return in_datab
		in_level = in_datab.description["Classification level"]
		if in_level not in ref.columns:
			print("Actual classification of database not contemplated in the file of reference")
			return in_datab
		if list(ref.columns).index(in_level) < list(ref.columns).index(out_level):
			print("Non si può ottenere una tassonomia più fine, si può solo salire nell'albero genealogico")
			return in_datab
		lil = {}
		for new in set(ref[out_level]):
			olds = ref[ref[out_level] == new][in_level]
			lil[new] = in_datab.samples[olds].sum(axis = 1)
		New = datab()	
		New.init_form_panda(pd.concat(lil, axis = 1), in_datab.description, in_datab.filename)
		return New
	
	
def Node_beet_centrality(data):	
	dik = nx.betweenness_centrality(data.graph)
	sort = sorted(list(dik), key = lambda d : dik[d])
	v = [dik[lab] for lab in sort]
	plt.figure(figsize=(9, 3))
	plt.subplot()
	plt.bar(sort, v)
	return [sort,v]

def Node_degree(data):
	dik = dict(data.graph.degree())
	sort = sorted(list(dik), key = lambda d : dik[d])
	v = [dik[lab] for lab in sort]
	plt.figure(figsize=(9, 3))
	plt.subplot()
	plt.bar(sort, v)
	return [sort,v]

def Edge_beet_centrality(data):
	dik = nx.edge_betweenness_centrality(data.graph)
	sort = sorted(list(dik), key = lambda d : dik[d])
	l = [str(lab[0]) + str(lab[1]) for lab in sort]
	v = [dik[lab] for lab in sort]
	plt.figure(figsize=(9, 3))
	plt.subplot()
	plt.bar(l, v)
	return [l,v]

def Jaccard_index_edges(data_1, data_2):
	inter = list(set(data_1.graph.edges()) & set(data_2.graph.edges()))
	nume = sum([min([data_1.graph.edges()[val]["weight"], data_2.graph.edges()[val]["weight"]]) for val in inter])
	deno = sum([data_1.graph.edges()[val]["weight"] for val in data_1.graph.edges()] + [data_2.graph.edges()[val]["weight"] for val in data_2.graph.edges()]) - nume
	return nume/deno
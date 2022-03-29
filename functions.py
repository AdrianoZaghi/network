import csv
import numpy
import math
from classes import datab, sample
import networkx as nx

#level can be "domain",	"phylum",	"class",	"order",	"family",	"genus" or	"otu"
def set_tax_level(in_datab, out_level, taxonomy_file):
	"""Reorganize the occurrence in different classifications level:
			in_datab 	 	datab object 	 	 a database containing the occurrency to be reorganized
			out_level 	 	string 	 	 	 	 a string thet have to match with the classification level you want as output
											     out_level level can be "domain",	"phylum",	"class",	"order",	"family",	"genus" or	"otu"
			taxonomy_file 	string 	 	 	 	 the filename of the CSV containing the informatino bout the correspondance between different classification levels
	Returns a different datab
	"""
	in_level = in_datab.description["Classification level"]
	new_samples = {}
	if out_level not in ["domain", "phylum", "class", "order", "family", "genus", "otu"]:
		print("nome livello sbagliato")
		return in_datab
	with open(taxonomy_file, mode = 'r', encoding = "utf-8") as r:
		ref = list(csv.DictReader(r))
		for s in in_datab.samples:
			new_occurrence = {}
			for raws in ref:
				Out, In = raws[out_level], raws[in_level]
				if Out in new_occurrence.keys():
					new_occurrence[Out] += in_datab.samples[s].occurrence[In]
				elif Out not in new_occurrence.keys():
					new_occurrence[Out] = in_datab.samples[s].occurrence[In]
			new_meta = in_datab.samples[s].metadata
			new_meta["Classification level"] = out_level
			new_samples[s] = sample(s, new_occurrence, new_meta)
	new_datab = datab()
	new_datab << [new_samples, in_datab.description.copy(), in_datab.filename]			
	new_datab.description["Classification level"] = out_level
	return new_datab	

#	New_tasso = set_tax_level(A, "phylum", "taxonomy.csv")
	
def filtering_prevalenza(data, preval):
	"""Filters input databases for the PREVALENZA
	An OTU (or other classification label) is discarded if has a prevalenza under second argument value (persentage)
	Returns a new datab
	"""
	in_labels = set(data.get_label_population())
	new_labels = {}
	for l in in_labels:
		new_labels[l] = 0
		for s in data.samples:
			if data.samples[s].occurrence[l] >= 1:
				new_labels[l] += 1
		if new_labels[l] < len(data.samples) * preval / 100:
			del new_labels[l]
	new_samples = {}
	for s in data.samples:
		new_occurrence = dict((k, data.samples[s].occurrence[k]) for k in new_labels)
		new_samples[s] = sample(s, new_occurrence, data.samples[s].metadata)
	new_data = datab()
	new_description = data.description.copy()
	new_description["Filtering procedure"] += ("Filtering for prevalenza at " + str(preval) +'%')
	new_data << [new_samples, new_description, data.filename]
	return new_data

#	pre = filtering_prevalenza(A, 10)

def filtering_mediana(data, mediana):
	"""Filters input databases fot the mediana
	An OTU (or other classification label) is discarded if has a mediana under second argument value (persentage)
	Returns a new datab	
	"""
	labels = set(data.get_label_population())
	new_labels = []
	for i in labels:
		top =  0
		bott = 6000
		for s in data.samples:
			if i in data.samples[s].occurrence.keys():
				if data.samples[s].occurrence[i] > top:
					top = data.samples[s].occurrence[i]
				if data.samples[s].occurrence[i] < bott:
					bott = data.samples[s].occurrence[i]
		if (top + bott) / 2 > mediana:
			new_labels.append(i)
	new_samples = {}
	for s in data.samples:		
		new_occurrence = dict((k, data.samples[s].occurrence[k]) for k in new_labels)
		new_sample = sample(s, new_occurrence, data.samples[s].metadata)
		new_samples[s] = new_sample
	new_description = data.description.copy()
	new_description["Filtering procedure"] += ("Filtering for mediana equal " + str(mediana))
	new_data = datab()
	new_data << [new_samples, new_description, data.filename]
	return new_data

#	med = filtering_mediana(A, 3)
#	filtrato = filtering_prevalenza(filtering_mediana(A, 5), 20)

def normalize(data, norm, normalization_name):
	"""This function normalize each sample of the datab with a given norm
		datab 	 	datab object 	 	 	the database to be normalized
		norm 	 	function name 	 	 	the function that will be used to normalize each sample
		normalization_name 	 	 	 	 	a string that will be added to the desctiption of the database
	If the database is already normalized, it's not modified and is returned as it is
	If not, is returned a new normalized database, without modifaing the existing one
	"""
	if data.description["Normalization"] == "nothing":
		new_samples = {}
		for s in data.samples:
			new_samples[s] = sample(s, norm(data.samples[s].occurrence), data.samples[s].metadata)
			new_samples[s].metadata["Normalization"] = normalization_name
		new_description = data.description.copy()
		new_description["Normalization"] = normalization_name
		new_data = datab()
		new_data << [new_samples, new_description, data.filename]
		return new_data
	else:
		print("Database già normalizzato")
		return data

#NORME

def L_1(s):
	"""A function that can be used to normalize sample occurrence with L1 norm
	"""
	new_occurrence = {}
	n = sum(s.values())
	for o in s:
		new_occurrence[o] = s[o] / n
	return new_occurrence

def CLR(s):
	"""A function that can be used to normalize sample occurrence with Center Log Rateo norm
	"""
	new_occurrence = {}
	P = 1
	for val in s:
		if s[val] > 0:
			P *= pow(s[val], 1/len(s))
	for o in s:
		if s[o] > 0:
			new_occurrence[o] = math.log(s[o] / P)
		else:
			new_occurrence[o] = 0
	return new_occurrence

# normato = normalize(A, L_1, "L_1")
	
def make_graph(c_matrix):
	"""Return a NetworkX object that represents a graph.
	It is obtained using as adjacency matrix the method c_matrix of the datab object in argument
	"""
	if len(c_matrix) == 0:
		print("La matrice è vuota")
		return
	return nx.from_numpy_array(c_matrix)

#	make_graph(A)

def edge_filtering(grap, density):
	"""From the NetworkX graph object in input, are removed the less hevi edges,
	untill the density sepcified is reached
	"""
	grap.remove_edges_from(nx.selfloop_edges(grap))
	sorted_by_weight = sorted(grap.edges(data="weight"), key=lambda tup: tup[2])
	while nx.density(grap) > density:
		grap.remove_edge(sorted_by_weight[0][0], sorted_by_weight[0][1])
		del sorted_by_weight[0]
	return grap
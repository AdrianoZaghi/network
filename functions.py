from classes import datab
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.lines as lines
import json
import numpy as np
import scipy.stats as st

#"domain","phylum","class","order","family","genus","otu"
def grow_tax_level(in_datab, out_level, taxonomy_file = "taxonomy.csv"):
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
	
def Json_for_cytoscape(data, filename):
	with open(filename + "_for_cytoscape.json", mode = "w", encoding ="utf-8") as out:
		cy = nx.cytoscape_data(data.graph)
		for edg in cy["elements"]["edges"]:
			edg["data"]["name"] = str(edg["data"]["source"])+" to "+str(edg["data"]["target"])
			edg["data"]["distance"] = abs(1/edg["data"]["weight"])
		b = nx.betweenness_centrality(data.graph)
		degree = dict(data.graph.degree())
		for node in cy["elements"]["nodes"]:
			node["data"]["bee centrality"] = b[node["data"]["name"]]
			node["data"]["degree"] = degree[node["data"]["name"]]
		print(json.dumps(cy), file = out)
	
def Otu_abundance(data):
	dik = dict(data.samples.sum(axis = 1))
	sort = sorted(list(dik), key = lambda d : dik[d])
	v = [dik[lab] for lab in sort]
	plt.figure(figsize=(9, 3))
	plt.subplot()
	plt.bar(sort, v)
	return [sort,v]

def get_filter_plot(data):
	sam = list(data.samples) 
	y = [data.samples.median()[l] for l in sam] 
	x = [list(data.samples[l]).count(0) for l in sam]
	colors = np.random.rand(3)
	somma = data.samples.sum(axis = 0)/50
	area = [somma[i] for i in sam]
	plt.scatter(x, y, s=area, c=colors, alpha=0.5)
	plt.yscale("log")
	plt.ylabel("Median")
	plt.xlabel("Number of 0 entries")
	plt.title("Filtering cuts")
	plt.plot([0,50],[5,5], color = "red")
	plt.plot([20,20],[0,1500], color = "red")
	plt.show()
	
def tre_test(A, B, C, mode = "degree"):
	if mode == "degree":
		a = [A.graph.degree[node] for node in A.graph.nodes]
		a = [0.5 if a[i] == 0 else a[i] for i in a]
		b = [B.graph.degree[node] for node in B.graph.nodes]
		b = [0.5 if b[i] == 0 else b[i] for i in b]
		c = [C.graph.degree[node] for node in C.graph.nodes]
		c = [0.5 if c[i] == 0 else c[i] for i in c]
#		print("AB " + str(st.ttest_ind(list(dict(A.graph.degree()).values()), list(dict(B.graph.degree()).values()))))
#		print("AC " + str(st.ttest_ind(list(dict(A.graph.degree()).values()), list(dict(C.graph.degree()).values()))))
#		print("BC " + str(st.ttest_ind(list(dict(B.graph.degree()).values()), list(dict(C.graph.degree()).values()))))
		print("AB " + str(st.chisquare(a, b, len(a)+len(b))))
		print("AC " + str(st.chisquare(a, c)))
		print("BC " + str(st.chisquare(b, c)))
#	elif mode == "centrality":
#		print("AB " + st.ttest_ind(dict(nx.betweenness_centrality(A.graph)).values(), dict(nx.betweenness_centrality(B.graph)).values()))
#		print("AC " + st.ttest_ind(dict(nx.betweenness_centrality(A.graph)).values(), dict(nx.betweenness_centrality(C.graph)).values()))
#		print("BC " + st.ttest_ind(dict(nx.betweenness_centrality(B.graph)).values(), dict(nx.betweenness_centrality(C.graph)).values()))
		
def histo_plot(data, mode = "degree"):
	fig, ax = plt.subplots()
	if mode == "degree":
		x = dict(data.graph.degree())
	elif mode =="centrality":
		x = nx.betweenness_centrality(data.graph)
	ax.hist([x[i] for i in x], 20)
	plt.xlabel(mode + " value")
	plt.ylabel("Occurrence")
	plt.title(mode + " distribution")
	
def comparison_network_features(data_a, data_b, mode = "centrality", X_lab = "L1", Y_lab = "SPARCC"):
	if mode == "degree":
		order = list(dict(data_a.graph.degree()))
		list_a = [data_a.graph.degree()[otu] for otu in order]
		list_b = [data_b.graph.degree()[otu] for otu in order]
	elif mode == "centrality":
		order = list(nx.betweenness_centrality(data_a.graph))
		list_a = [nx.betweenness_centrality(data_a.graph)[otu] for otu in order]
		list_b = [nx.betweenness_centrality(data_b.graph)[otu] for otu in order]
	elif mode == "correlation":
		order = list(data_a.graph.edges)
		list_a = [data_a.graph.edges[otus]["weight"] for otus in order]
		list_b = [data_b.graph.edges[otus]["weight"] for otus in order]
	plt.xlabel(X_lab + " " + mode)
	plt.ylabel(Y_lab + " " + mode)
	plt.plot(list_a, list_b, 'o', color='purple')
	plt.plot([min(list_a),max(list_a)],[min(list_b),max(list_b)], color = "black")
	plt.show()
	
	
def make_tabella(dict_L1, dict_CLR, dict_SP):
	sort_L1 = sorted(list(dict_L1), key = lambda d : abs(dict_L1[d]["weight"]), reverse=True)[0:10]
	sort_CLR = 	sorted(list(dict_CLR), key = lambda d : abs(dict_CLR[d]["weight"]), reverse=True)[0:10]
	sort_SP = sorted(list(dict_SP), key = lambda d : abs(dict_SP[d]["weight"]), reverse=True)[0:10]
	righe = set(sort_CLR + sort_L1 + sort_SP)
	with open("latex_tabel.tex", mode = 'w', encoding = "utf-8") as f:
		print("\\begin{center} \n \\begin{tabular}{ |c|c|c|c| }\n", file = f)
		print("  & L1 & CLR & SPARCC\\\\\n\\hline", file = f)
		for r in righe:
			print(r[0]+" to "+r[1], end = " & ", file = f)
			if r in sort_L1:
				print(round(dict_L1[r]["weight"], 3), end = " & ", file = f)
			else:
				print(" - ", end = " & ", file = f)
			if r in sort_CLR:
				print(round(dict_CLR[r]["weight"], 3), end = " & ", file = f)
			else:
				print(" - ", end = " & ", file = f)
			if r in sort_SP:
				print(round(dict_SP[r]["weight"], 3), end = " \\\\\n", file = f)
			else:
				print(" - ", end = " \\\\\n", file = f)
		print("\\end{tabular}\n\\end{center}", file = f)
				
def last(A, B, step, metod_1, metod_2):
	edges_a = list(A.graph.edges)
	edges_a = sorted(edges_a, key = lambda lab : abs(A.graph.edges[lab]["weight"]))
	edges_b = list(B.graph.edges)
	edges_b = sorted(edges_b, key = lambda lab : abs(B.graph.edges[lab]["weight"]))
	if len(edges_a) != len(edges_b) or len(edges_a) != (len(A.graph.nodes) * (len(A.graph.nodes)-1))/2 :
		print("qualcosa è andato storto")
		return
	color = []
	for i in range(0, len(edges_a), step):
		a = edges_a[i:]
		raw = []
		for j in range(0, len(edges_b), step):
			b = edges_b[j:]
			b_inv = [(edge[1], edge[0]) for edge in b]
			common = len(list(set(a).intersection(b + b_inv)))
			raw.append((common*common) / (len(a)*len(b)) - common / max(len(a), len(b)))
#			raw.append((common*common) / (len(a)*len(b)))
		color.append(raw)
	plt.imshow(color, interpolation='bilinear')
	plt.xlabel("Links removed from " + metod_1)
	plt.ylabel("Links removed from " + metod_2)
	plt.colorbar()
	plt.show()
	return color
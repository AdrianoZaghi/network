import os
import pandas as pd
import numpy as np
import networkx as nx

class datab:
	def __init__(self, file_name = None):
		self.c_matrix = None
		self.graph = None
		if file_name == None:
			self.description = {}
			self.filename = ""
			self.samples = None
		else:
			self.samples = pd.read_csv(file_name, index_col=0)
			self.filename = file_name.split('.')[0]
			self.description = {"Classification level":"otu", "Filtering procedure" : "", "Normalization": "nothing", "Square?" : False}
	
	def __lshift__(self, dataframe_description_filename):	#	<< operator
		self.samples = dataframe_description_filename[0]
		self.description = dataframe_description_filename[1]
		self.filename = dataframe_description_filename[2]
	
	def get_info(self):
		"""Print on the terminal some usefull informations about the database"""
		print(self.samples.info())
		print(self.description)
		
	def filter_median(self, m):
		not_keep = dict(self.samples.median() < m)
		k = [lab for lab in not_keep if not_keep[lab]]
		self.samples = self.samples.drop(k, axis = "columns")
		
	def filter_prevalence(self, p):
		limit = (1-p)*self.samples.shape[0]
		not_keep = [l for l in self.samples if list(self.samples[l]).count(0) > limit]
		self.samples = self.samples.drop(not_keep, axis = "columns")
	
	def get_sparcc_matrix(self, iterations):
		if self.description["Normalization"] != "nothing":
			print("Per questo metodo, i dati non devono essere normalizzati")
			return
		name = self.filename + "_sparcc_utility.tsv"
		with open(name, mode = 'w', encoding="utf-8") as tab:
			self.samples.T.to_csv(tab, sep = '\t')
		os.system("python3 ../SparCC3/SparCC.py " + name + " -i " + str(iterations) + " --cor_file=" + name.split("_utility")[0] + "_cor_matrix.tsv")
		self.c_matrix = pd.read_csv(name.split("_utility")[0] + "_cor_matrix.tsv", sep='\t', index_col=0).to_numpy()
		return self.c_matrix
		
	def get_pearson_matrix(self, mode = "L1"):
		if mode == "L1":
			norm = self.samples.sum(axis = 1)
			self.samples.div(norm, axis = 0)
		elif mode == "CLR":
			temporary = self.samples.replace(0, np.nan)
			G = pow(temporary.product(axis = 1, skipna = True), 1/len(self.sampls.shape[1]))
			self.samples = np.log(self.samples.div(G, axis = 0))
		else:
			print("normalization mode not implemented")
			return 0
		self.description["Normalization"] = mode
		self.c_matrix = np.corrcoef(np.array(self.samples), rowvar = False)
		return self.c_matrix
	
	def make_graph(self, density):
		self.graph = nx.from_numpy_array(self.c_matrix)
		self.graph.remove_edges_from(nx.selfloop_edges(self.graph))
		sorted_by_weight = sorted(abs(self.graph.edges(data="weight")), key = lambda tup: abs(tup[2]))
		while nx.density(self.graph) > density:
			self.graph.remove_edge(sorted_by_weight[0][0], sorted_by_weight[0][1])
			del sorted_by_weight[0]
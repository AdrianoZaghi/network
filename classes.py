import os
import pandas as pd
import numpy as np
import networkx as nx

class datab:
	"""Datab class is made to contain the data from a csv file containing the otu (columns) occurrence of metagenomic samoles (raws)
	It's methods are:
		samples 	 	 	 	 	 A pandas dataframe containing the data readed forom the csv file
		filename 	 	 	 	 	 A string that keep track of the file name used to construct the object
		description 	 	 	 	 A dictionary with various labeled informations about the operations performed on the object
			"Classification level" 	 	Specify if the occurrence is orderedo byo otu or by an other classification (phylim, family, class, order...)+
			"Filtering procedure" 	 	How the database have been filtered
			"Normalization" 	 	 	How the databas have been normalized
			"Correlation" 	 	 	 	How the correlation matrix have been obrained
			"Graf filtering procedure" 	How the graph have been filtered
		c_matrix 	 	 	 	 	 A pandas matrix, containing the correlation matrix evaluated from the data
		graph 	 	 	 	 	 	 A networkx graph that is constructed from the correlation matrix
	"""
	def __init__(self, file_name = None):
		"""The constructor can have at most 1 argument. If not provided, the database is created empthy object
			file_name 	 	 path to the occurrence file
		"""
		self.c_matrix = None
		self.graph = None
		if file_name == None:
			self.description = {}
			self.filename = ""
			self.samples = None
		else:
			self.samples = pd.read_csv(file_name, index_col=0)
			self.filename = file_name.split('.')[0]
			self.description = {"Classification level":"otu", "Filtering procedure" : "", "Normalization": "nothing", "Correlation": "not yet", "Grph filtering procedure": "not yet"}
	
	def init_form_panda(self, samples, description, filename = "no_file"):
		"""Initialize empthy datab objects with from runtime variables
			samples 	 	a pandas dataframe containing occurrence
			description 	a description for the new datab
			filename 	 	a filename can be rovided, if is intended to refer the data to an external database
		"""
		self.samples = samples
		self.description = description
		self.filename = filename
	
	def get_info(self):
		"""Print on the terminal some usefull informations about the database
		"""
		graph_info = {}
		print("Samples informations")
		print(self.samples.info())
		print("Pipline history")
		print(self.description)
		if self.graph == None:
			print("Grafico non ancora costruito")
		else:
			print("Returno le seguenti informazioni sul network")
			print("Nodes info : {names : betweenness centrality}")
			graph_info["Nodes"] = nx.beetweennes_centrality(self.graph)
			graph_info["Edges"] = nx.edge_betweenness_centrality(self.graph)
		return graph_info
				
	def filter_median(self, m):
		"""Remove all the database columns which have a median value under the value of the argument m
		"""
		not_keep = dict(self.samples.median() < m)
		k = [lab for lab in not_keep if not_keep[lab]]
		r = self.samples[k].copy()
		self.samples = self.samples.drop(k, axis = "columns")
		self.description["Filtering procedure"] += "Filtered out clumn with median under " + str(m) + "\n" 
		return r
		
	def filter_prevalence(self, p):
		"""Remove all the database columns which have less the p% of non 0 values
		"""
		limit = (100-p)*self.samples.shape[0]/100
		not_keep = [l for l in self.samples if list(self.samples[l]).count(0) > limit]
		self.samples = self.samples.drop(not_keep, axis = "columns")
		self.description["Filtering procedure"] += "Filtered out column with more then " + str(p) + "% 0 values\n"
	
	def get_sparcc_matrix(self, iterations):
		"""Initialize the method c_matrix with a correlation matrix obtained with sparcc algorithm
		The algorithm performs "iterations" iterations
		The function also returns the matrix
		"""
		if self.description["Normalization"] != "nothing":
			print("Per questo metodo, i dati non devono essere normalizzati")
			return
		name = self.filename + "_sparcc_utility.tsv"
		with open(name, mode = 'w', encoding="utf-8") as tab:
			self.samples.T.to_csv(tab, sep = '\t')
		os.system("python3 ../SparCC3/SparCC.py " + name + " -i " + str(iterations) + " --cor_file=" + name.split("_utility")[0] + "_cor_matrix.tsv")
		self.c_matrix = pd.read_csv(name.split("_utility")[0] + "_cor_matrix.tsv", sep='\t', index_col=0)
		self.description["Correlation"] = "Correelation evaluated with SparCC"
		return self.c_matrix
		
	def get_pearson_matrix(self, mode = "L1"):
		"""Initialize the method c_matrix with a Pearson correlation matrix
		This function ofers 2 different modes to normalize the datas before computing correlation:
			mode = "L1" (default) 	 	Normalize the data of each sample for their sum
			mode = "CLR" 	 	 	 	Normalize the data of each sample for their center log rateo
		The function also returns the matrix
		"""
		if mode == "L1":
			norm = self.samples.sum(axis = 1)
			self.samples = self.samples.div(norm, axis = 0)
			self.description["Normalization"] = "Normalized with L1"
		elif mode == "CLR":
	#		temporary = self.samples.replace(0, np.nan)
			temporary = self.samples.replace(to_replace = 0, value = 0.5)
			G = pow(temporary.product(axis = 1, skipna = True), 1/self.samples.shape[1])
			self.samples = np.log(self.samples.div(G, axis = 0))
			self.description["Normalization"] = "Normalized with CLR"
		else:
			print("normalization mode not implemented")
			return 0
		self.description["Normalization"] = mode
		self.c_matrix = self.samples.corr(method = "pearson")
		self.description["Correlation"] = "Correlation evalueted with Pearson"
		return self.c_matrix
	
	def make_graph(self, density):
		"""Initialize the graph method with a networkx object obtained using c_matrix as an adjacency matrix
		The graph is weighted and undirected
		After the graph is created, self loops are removed
		Then are also removed the smallest edges untill is reached the density specified in the argument
		"""
		self.graph = nx.from_pandas_adjacency(self.c_matrix)
		self.graph.remove_edges_from(nx.selfloop_edges(self.graph))
		sorted_by_weight = sorted(list(self.graph.edges(data="weight")), key = lambda tup: abs(tup[2]))
		while nx.density(self.graph) > density:
			self.graph.remove_edge(sorted_by_weight[0][0], sorted_by_weight[0][1])
			del sorted_by_weight[0]
		self.description["Graph filtering procedure"] = "Removing self loops\nRemoving edges untill i get the density of " + str(density) + "\n"

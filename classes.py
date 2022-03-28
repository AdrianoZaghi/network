import csv
import os
import pandas as pd
import numpy as np

class sample:
	"""Class used to collect sample informations:
			name 	 	string 	 	 	should be the name of the sample
			occurrence 	dicctionary 	should collect the occurrency of all the different kind of OTU (or other classification), labeled with their name
			metadata 	dictionary 	 	should collect the additional informations about the sample contained in the metadata file
	"""
	def __init__(self, name, occurrence, metadata):
		self.name = name
		self.occurrence = occurrence.copy()
		self.metadata = metadata.copy()

class datab:
	"""Class used as a container of multiple samples, should represent the database used
			samples 	 	dictionary 	 should contain all the samples, labeled with their name
			description 	dictionary 	 contains various fields that are modified when the database is manipulated. Keeps track of all the poerations made on it
				"Classification level" 	: 	The occurrence are expressed intherms of OTU, phyum, genus, ....
				"Filtering procedure" 	: 	What are the chriterion used to filter the OTUs (or other level) before performing any operation
				"Normalization" 	 	: 	What mathod is uded to normalize the database ebtries
				"Square" 	 	 	 	: 	True if all the samples have the same se of labels
			c_matrix 	 	 2D ndarray  contains the value of correlation between the various kind of bacteria that are evaluted from all the database
			filename 	 	 string 	 keep track of the filename the occurrencies are taken from"""
	def __init__(self, file_name = None, metadata_file_name = None):
		"""__init__() 	 	 it can recive 2 arguments or none:
			file_name 	 	 	 should be a csv file with the occurrence of different OTUs
			metadata_file_name   should be a csv file conitaining all other infos about the samples"""
		self.c_matrix = []
		self.samples = {}
		if file_name == None:
			self.description = {}
			self.filename = ""
		else:
			self.filename = file_name.split('.')[0]
			self.description = {"Classification level":"otu", "Filtering procedure" : "", "Normalization": "nothing", "Square?" : False}
			with open(file_name, encoding= "utf-8", mode = "r") as d, open(metadata_file_name, encoding="utf-8", mode ="r") as e:
				csv_dati = csv.DictReader(d)
				csv_metadati = csv.DictReader(e)
				for righe in csv_dati:
					occ = righe
					name = occ.pop("")
					for elem in occ:
						occ[elem] = int(occ[elem])
					for righe in csv_metadati:
						if righe[""] == name:
							met = righe
							break
					del met[""]
					met["Classification level"] = "otu"
					self.samples[name] = sample(name, occ, met)

	def __lshift__(self, list_of_2_dict_and_the_name):	#	<< operator
		self.samples, self.description, self.filename = list_of_2_dict_and_the_name

	def encert_square(self):
		first_sample_len = len(list(self.samples.values())[0].occurrence)
		r = all(len(self.samples[s].occurrence) == first_sample_len for s in self.samples)
		self.description["Square?"] = r
		return r

	def get_label_population(self):
		labels = []
		for s in self.samples:
	 		labels = labels + list(self.samples[s].occurrence.keys())
		return labels
	
	def get_info(self):
		"""Print on the terminal some usefull informations about the database"""
		print(self.description)
		print("ci sono ", len(self.samples), " campioni")
		print("elenco i campioni con annessa quantità di entrate")
		for s in self.samples:
			print(s, "      ", len(self.samples[s].occurrence))
			
	def print_tsv_for_sparcc(self, file_name):
		with open(file_name, mode = 'w', encoding="utf-8") as tab:
			print("", end='', file = tab)
			label_list = list(list(self.samples.values())[0].occurrence.keys())
			for s in self.samples:
				print("\t\""+s+"\"", end='', file = tab)
			print("", file = tab)
			for l in label_list:
				print("\""+l+"\"", end='', file = tab)
				for s in self.samples:
					print("\t", self.samples[s].occurrence[l], end='', file=tab)
				print("", file= tab)
	
	def get_sparcc_matrix(self, iterations):
		"""If the data haven't be normalized anyhouw and is square, this function initialize the method c_matrix with the correlation values obtained from sparcc algorithm
		Moreover it retruns this matrix
		The argument is the number of iterations that the algorithm is going to perform
		If the conditions are not observed, the function prints an error message and retuns nothing without initializing enithing
		"""
		self.encert_square()
		if self.description["Square?"] == False:
			print("Occorrenze non quadrate")
			return
		if self.description["Normalization"] != "nothing":
			print("Per questo metodo, i dati non devono essere normalizzati")
			return
		name = self.filename + "_sparcc_utility.tsv"
		self.print_tsv_for_sparcc(name)
		os.system("python3 ../SparCC3/SparCC.py " + name + " -i " + str(iterations) + " --cor_file=" + name.split('.')[0] + "_cor_matrix.tsv")
		self.c_matrix = pd.read_csv(name.split('.')[0] + "_cor_matrix.tsv", sep='\t')
		return self.c_matrix
		
	def get_pearson_matrix(self):
		"""This function evluate the c_matrix with pearson correlation and retunrs it.
		The data have to be squared, otherwise it prints an error message and don't do anything else
		Before doing the caluclation, it prints out the way in which the databade have been normalized
			The data should always be normalized before evaluating pearson correlation
		"""
		self.encert_square()
		if self.description["Square?"] == False:
			print("Occorrenze non quadrate")
			return
		print("The dataset has been previously normalized with " + self.description["Normalization"])
		label_list = list(set(self.get_label_population()))
		sample_ids = list(self.samples.keys())
		stupid = []
		for s in sample_ids:
			List = []
			for o in label_list:
				List.append(self.samples[s].occurrence[o])
			stupid.append(List)
		self.c_matrix = np.corrcoef(np.array(stupid), rowvar = False)
		return self.c_matrix
		
	
#	A = datab("abundance.csv","config.json","metadata.csv")
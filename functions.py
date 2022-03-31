from classes import datab
import pandas as pd


def grow_tax_level(in_datab, out_level, taxonomy_file):
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
		New = pd.concat(lil, axis = 1)
		return New
			
A = datab("abundance.csv")
grow_tax_level(A, "order", "taxonomy.csv")
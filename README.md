# README

classes.py and functions.py provides tools to work with csv files containing samples of bacteira popultion.
In each raw should be contained the occurrence of all the different kind of bacretia that appeared in the sample.
Bacteria devision is suposed to be in therms of otu.
Classes and functions are thougt to be used on the ipython console, and various kind of operation are supported:

- Initializing a datasbase with the data in the "occorrenze.csv"\
	A = datab("occorrence.csv")
	
- Changing classification level to an higher order classification(from otu to phylum for example)\
A taxonomy file that specifies the correspndances between different levelso fo classification has to be provided (see taxonomy.csv as an example)
	A_1 = set_tax_level(A, out_level, "taxonomy.csv")\
	out_level can be: "domain",	"phylum",	"class",	"order",	"family",	"genus" or	"otu"
	
- Filtering the database and removing low occurrency otus\
	- Filtering out all the otu that have prevalenza under cetrain LEVEL in the database (percentage)\
	A.filtering_prevalenza(LEVEL)
	
	- Filtering out all the otu that have mediana under cetrain LEVEL in the database\
	A.filtering_mediana(LEVEL)
	
- Creating a correlation matrix of samples occurrency in various ways:\
	- Pearson correlation with "L1" data normalization\
	A.get_pearson_matrix("L1")
	
	- Pearson correlation with "CLR" data normalization\
	A.get_pearson_correlation("CLR")
	
	- SparCC correlation matrix. In argument can be speciied the number of iterations of such algorithm\
	A.get_sparcc_matrix(iterations)
	
- Making a graph (NetworkX object) from this correlation matrix. In the construction are removed self loops and also the less hevy links, untill the specified link density is obtained\
	A.make_graph(density)

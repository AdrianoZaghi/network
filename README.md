# README

classes.py and functions.py provides tools to work with csv files containing samples of bacteira popultion.
In each raw should be contained the occurrence of all the different kind of bacretia that appeared in the sample
Bacteria devision is suposed to be in therms of otu
Should be provided also an additional metadata.csv file, that can be mepthy, but can also contain various kind of informations about the samples
Classes and functions are thougt to be used on the ipython console, and various kind of operation are supported:

- Initializing a datasbase with the data in the "occorrenze.csv" file and those in the "metadata.csv" file
	A = datab("occorrence.csv", "metadata.csv")
	
- Changing classification level (from otu to phylum for example)
	A_1 = set_tax_level(A, out_level, "taxonomy.csv")
	out_level can be: "domain",	"phylum",	"class",	"order",	"family",	"genus" or	"otu"
	
- Filtering the database and removing low occurrency otus
	- Filtering out all the otu that have prevalenza under cetrain LEVEL in the database (percentage)
	A_1 = filtering_prevalenza(A, LEVEL)
	
	- Filtering out all the otu that have mediana under cetrain LEVEL in the database
	A_1 = filtering_mediana(A, LEVEL)
	
- Normalize the samples in the database in various ways
	- L_1 normalization
	A_1 = normalize(A, L_1, "L1")
	
	- Center Log Rateo (CRL) normalization
	A_1 = normalize(A, CLR, "Center Log Rateo")
	
- Evaluate the correlation network that can be obtained from the database, between otu occcurrency
	- Pearson correlation matrix **from previously normalized database**
	make_graph(A)
	
	- SparCC correlation matrix with **from previously non normalized database**. SparCC algorithm will iterate *n* times
	A.get_sparcc_matrix(*n*)
	make_graph(A)
	
- Remove link less hevy links form the correlation network untill a desired luink *density* is obtained
	- edge_filtering(A, *density*)

import classes as c
import functions as f

		
#creazione del grafo in forma di network
A = c.datab("abundance.csv")
B = c.datab("abundance.csv")
A.filter_prevalence(20)
B.filter_prevalence(20)
A.filter_median(5)
B.filter_median(5)
#...
A.get_sparcc_matrix(20) 		 
B.get_pearson_matrix("L1") 	 	# 	A.get_pearson_matrix("CLR")
#...
A.make_graph(0.05)
B.make_graph(0.05)
#
uno = f.Node_beet_centrality(A)
due = f.Edge_beet_centrality(A)
#
J = f.Jaccard_index_edges(A, B)
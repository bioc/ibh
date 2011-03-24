\name{ibh-package}
\Rdversion{1.1}
\alias{ibh-package}
\docType{package}
\title{
Interaction Based Homogeneity
}
\description{
This package contains methods for evaluating Interaction Based Homogeneity for lists of genes. Given a gene list of n genes, 
we first form an adjacency matrix A whose rows and columns are genes in the list where \deqn{A_{ij} = 1} if genes i and j 
have an interaction in the network and \deqn{A_{ij}=0} otherwise. The Interaction Based Homogeneity for a gene list 
\deqn{L=\{g_1, g_2, ..., g_n\}} of size n is then calculated as:

\deqn{InteractionBasedHomogeneity(L) ={{\sum_{i=1}^n\sum_{j=1}^n{A_{ij}}}\over{n^2}}}
}
\details{
\tabular{ll}{
Package: \tab ibh\cr
Type: \tab Package\cr
Version: \tab 1.0.0\cr
Date: \tab 2011-01-19\cr
License: \tab GPL (version 2 or newer)\cr
LazyLoad: \tab yes\cr
}
The user can provide his own interaction list or can use predefined gene lists which are created based on the BioGRID Interactions. 
Both gene lists and result of clustering methods such as kmeans or hclust can be used as inputs. Entrez identifiers, unique ids
(systematic names) or official names can be used as gene/protein identifiers.
}
\author{
Kircicegi Korkmaz <e102771@ceng.metu.edu.tr>, Volkan Atalay <volkan@ceng.metu.edu.tr>, Rengul Cetin-Atalay <rengul@bilkent.edu.tr>

Maintainer: Kircicegi Korkmaz <e102771@ceng.metu.edu.tr>}
\references{
Stark C, Breitkreutz BJ, Reguly T, Boucher L, Breitkreutz A, 
Tyers M. \emph{Biogrid: A General Repository for Interaction Datasets.} Nucleic Acids Res. Jan1; 34:D535-9}
\keyword{ package }
\examples{
	require(simpIntLists)
	data(ArabidopsisBioGRIDInteractionEntrezId)
	listofGeneList <- list(list(839226,817241, 824340, 832179, 818561, 831145, 838782, 826404), 
					list( 832018, 839226, 839226, 838824));
	ibhForMultipleGeneLists(ArabidopsisBioGRIDInteractionEntrezId, 
			listofGeneList)

	require(simpIntLists)
	listofGeneList <- list(list(839226,817241, 824340, 832179, 818561, 831145, 838782, 826404), 
list( 832018, 
				839226, 839226, 838824));
	ibhForMultipleGeneListsBioGRID(listofGeneList, 
				 organism="arabidopsis",
				 idType = "EntrezId");
	listofGeneList <- list(list("YJR151C", "YBL032W",   "YAL040C",   "YBL072C",   "YCL050C",
					   "YCR009C"), 
				list("YDR063W", "YDR074W", "YDR080W", "YDR247W", "YGR183C", "YHL033C"), 
				list("YOL068C" , "YOL015W" ,  "YOL009C" ,  "YOL004W" ,  "YOR065W"  ));
	ibhForMultipleGeneListsBioGRID(listofGeneList, organism="yeast", 
			idType = "UniqueId");


	require(yeastCC)
	require(stats)
	data(yeastCC)
	subset <- exprs(yeastCC)[1:50,]
	d <- dist(subset,method="euclidean")
	k <- kmeans(d, 3);
	ibhClusterEvalBioGRID(k$cluster, rownames(subset),
			organism="yeast", idType="UniqueId")

}

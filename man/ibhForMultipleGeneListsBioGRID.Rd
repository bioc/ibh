\name{ibhForMultipleGeneListsBioGRID}
\Rdversion{1.1}
\alias{ibhForMultipleGeneListsBioGRID}
\title{
Calculate interaction based homogeneity for multiple gene lists according to the BioGRID Interactions}
\description{
This function calculated interactios based homogeneity for multiple gene lists according to the  BioGRID 
Interactions for seven organisms: Arabidopsis thaliana(arabidopsis), Caenerhabditis elegans(c.elegans), 
Drosophila melanogaster (fruitFly), Homo sapiens (human), Mus musculus (mouse), Saccharomyces cerevisae 
(yeast), Schizosaccharomyces pombe (s.pombe). Unique ids(systematic names), official names or Entrez ids 
can be used as identifier type.}
\usage{
ibhForMultipleGeneListsBioGRID(listofGeneList,
			 organism, idType = "EntrezId")
}
\arguments{
  \item{listofGeneList}{
List of list of genes/proteins for which interaction based homogeneity is evaluated.}
  \item{organism}{
  Organism name. Can be one of 'arabidopsis', 'c.elegans', 'fruitFly', 'human', 'mouse', 'yeast', 's.pombe'. 
}
  \item{idType}{
Type of identifier. Can be one of 'EntrezId', 'Official' and 'UniqueId'.
}
}
\value{
 A vector of floats representing interaction based homogeneity for each list
}
\references{
Stark C, Breitkreutz BJ, Reguly T, Boucher L, Breitkreutz A, 
Tyers M. \emph{Biogrid: A General Repository for Interaction Datasets.} Nucleic Acids Res. Jan1; 34:D535-9}

 \seealso{
 \code{\link{ibh}}}
\examples{
	require(simpIntLists)
	listofGeneList <- list(list(839226,817241, 824340, 832179, 818561, 831145, 
				838782, 826404), 
				list( 832018, 839226, 839226, 838824));
	ibhForMultipleGeneListsBioGRID(listofGeneList, 
			organism="arabidopsis",  idType = "EntrezId");
	listofGeneList <- list(list("YJR151C", "YBL032W",   "YAL040C",   "YBL072C",   
				"YCL050C",   
				"YCR009C"), list("YDR063W", "YDR074W", "YDR080W", "YDR247W", 
				"YGR183C", "YHL033C"), list("YOL068C" , "YOL015W" ,  
				"YOL009C" , "YOL004W" , "YOR065W"  ));
	ibhForMultipleGeneListsBioGRID(listofGeneList, 
		organism="yeast", idType = "UniqueId");
}
\keyword{ file}

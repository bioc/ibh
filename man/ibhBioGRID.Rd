\name{ibhBioGRID}
\Rdversion{1.1}
\alias{ibhBioGRID}
\title{
Calculate interaction based homogeneity for a gene list according to the BioGRID Interactions}
\description{
This function calculated interactios based homogeneity for a gene list according to the  BioGRID 
Interactions for seven organisms: Arabidopsis thaliana(arabidopsis), Caenerhabditis elegans(c.elegans),
Drosophila melanogaster (fruitFly), Homo sapiens (human), Mus musculus (mouse), Saccharomyces cerevisae 
(yeast), Schizosaccharomyces pombe (s.pombe). Unique ids(systematic names), official names or Entrez ids 
can be used as identifier type.
}
\usage{
ibhBioGRID(geneList, organism, idType = "EntrezId")
}
\arguments{
  \item{geneList}{
list of genes/proteins for which interaction based homogeneity is evaluated}
  \item{organism}{
  Organism name. Can be one of 'arabidopsis', 'c.elegans', 'fruitFly', 'human', 'mouse', 'yeast', 's.pombe'. 
}
  \item{idType}{
Type of identifier used. Can be one of 'EntrezId', 'Official' and 'UniqueId'}
}
\value{
 Float representing interaction based homogeneity for each list
}
\references{
Stark C, Breitkreutz BJ, Reguly T, Boucher L, Breitkreutz A, 
Tyers M. \emph{Biogrid: A General Repository for Interaction Datasets.} Nucleic Acids Res. Jan1; 34:D535-9}
\seealso{
\code{\link{ibh}},\code{\link{ibhForMultipleGeneListsBioGRID}}
}
\examples{
	require(simpIntLists)
	geneList <-  list(839226,817241, 824340, 832179, 818561, 831145, 838782, 826404);
	ibhBioGRID(geneList, organism="arabidopsis",
						idType = "EntrezId");
	geneList <- list("YJR151C", "YBL032W",   "YAL040C",   "YBL072C",   "YCL050C",   "YCR009C");
	ibhBioGRID(geneList, organism="yeast", idType = "UniqueId");
}
\keyword{ file}
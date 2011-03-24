\name{ibhClusterEvalBioGRID}
\Rdversion{1.1}
\alias{ibhClusterEvalBioGRID}
\title{
Evaluate clustering with interaction based homogeneity using BioGRID Interactions
}
\description{
This function calculated interaction based homogeneity for a clustering result to the  BioGRID 
Interactions for seven organisms: Arabidopsis thaliana(arabidopsis), Caenerhabditis elegans(c.elegans), 
Drosophila melanogaster (fruitFly), Homo sapiens (human), Mus musculus (mouse), Saccharomyces cerevisae 
(yeast), Schizosaccharomyces pombe (s.pombe). Unique ids(systematic names), official names or Entrez ids 
can be used as identifier type.}
\usage{
ibhClusterEvalBioGRID(cluster, allGenesList, 
    organism, idType = "EntrezId")
}
\arguments{
  \item{cluster}{
result of clustering
}
  \item{allGenesList}{
list of genes in the same order of clustering object
}
  \item{organism}{
organism name, can be one of 'arabidopsis', 'c.elegans', 'fruitFly', 'human', 'mouse', 'yeast', 's.pombe'. 
}
  \item{idType}{
type of identifier, can be one of 'EntrezId', 'Official' and 'UniqueId'.
}
}
\value{
 A vector of floats representing interaction based homogeneity for each cluste
}
\references{
Stark C, Breitkreutz BJ, Reguly T, Boucher L, Breitkreutz A, 
Tyers M. \emph{Biogrid: A General Repository for Interaction Datasets.} Nucleic Acids Res. Jan1; 34:D535-9}

\examples{
require(yeastCC)
require(stats)
require(simpIntLists)
data(yeastCC)
subset <- exprs(yeastCC)[1:50,]
d <- dist(subset,method="euclidean")
k <- kmeans(d, 3);
ibhClusterEvalBioGRID(k$cluster, rownames(subset), 
			organism="yeast", idType="UniqueId")
}
\keyword{ file}
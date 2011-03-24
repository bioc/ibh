\name{ibhClusterEval}
\Rdversion{1.1}
\alias{ibhClusterEval}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Evaluate clustering with interaction based homogeneity}
\description{
This function calculated interaction based homogeneity for a clustering result. }
\usage{
ibhClusterEval(cluster, allGenesList, interactionList)
}
\arguments{
  \item{cluster}{
result of clustering
}
  \item{allGenesList}{
list of genes in the same order of clustering object
}
  \item{interactionList}{
list containing the interactions. For each gene/protein, the is an entry in the list with "name" 
containing name of the gen/protein and "interactors" containing the list of genes/proteins interacting with it.}
}
\value{
 A vector of floats representing interaction based homogeneity for each cluste
}
\examples{
require(yeastCC)
require(stats)
data(yeastCC)
require(simpIntLists)
data(YeastBioGRIDInteractionUniqueId)

subset <- exprs(yeastCC)[1:50,]
d <- dist(subset,method="euclidean")
k <- kmeans(d, 3);
ibhClusterEval(k$cluster, rownames(subset), 
			YeastBioGRIDInteractionUniqueId)
  }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ file}

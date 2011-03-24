\name{ibh}
\alias{ibh}
\title{ Calculate interaction based homogeneity for the given gene list according to the interaction list}
\description{
  This function calculates interaction based homogeneity of the given gene list according to the 
interaction list
}
\usage{
ibh(interactionList, geneList)}
\arguments{  
  \item{interactionList}{List containing the interactions. For each gene/protein, the is an entry in the list 
with "name" containing name of the gen/protein and "interactors" containing the list of genes/proteins 
interacting with it.}
  \item{geneList}{List of genes/proteins for which interaction based homogeneity is evaluated.}
}
\value{ Interaction based homogeneity value as float
}
\author{Kircicegi Korkmaz}
\seealso{  \code{\link{ibhForMultipleGeneLists}}}
\examples{
	require(simpIntLists)
	data(ArabidopsisBioGRIDInteractionEntrezId)
	geneList <-  list(839226,817241, 824340, 832179, 818561, 831145, 838782, 826404);
	ibh(ArabidopsisBioGRIDInteractionEntrezId,geneList);

}
\keyword{file}
\name{ibhForMultipleGeneLists}
\alias{ibhForMultipleGeneLists}
\title{ Calculate interaction based homogeneity for multiple gene lists according to the interaction list}
\description{
  This function calculated interactios based homogeneity for multiple gene lists according to the interaction list
}
\usage{
ibhForMultipleGeneLists(interactionList, listofGeneList)}
\arguments{  
  \item{interactionList}{List containing the interactions. For each gene/protein, the is an entry in the list with 
"name" containing name of the gen/protein and "interactors" containing the list of genes/proteins interacting with it.}
  \item{listofGeneList}{List of list of genes/proteins for which interaction based homogeneity is evaluated. }
}
\value{ A vector of floats representing interaction based homogeneity for each list
}
\author{Kircicegi Korkmaz}
\seealso{  \code{\link{ibh}}}
\examples{
	require(simpIntLists)
	data(ArabidopsisBioGRIDInteractionEntrezId)
	listofGeneList <- list(list(839226,817241, 824340, 832179, 818561, 
					831145, 838782, 826404), 
                        list( 832018, 839226, 839226, 838824));
	ibhForMultipleGeneLists(
			ArabidopsisBioGRIDInteractionEntrezId, listofGeneList)

}
\keyword{file}
% --- Source file: readUndirectedInteractionsFromCsv.Rd ---
\name{readUndirectedInteractionsFromCsv}
\alias{readUndirectedInteractionsFromCsv}
\title{ Read undirected interactions from csv and create the interaction list}
\description{
  This function reads the undirected interactions from a csv file and creates the interaction list. The csv file
 must contain two names: first gene/protein name, second the interactor.
}
\usage{
readUndirectedInteractionsFromCsv(fileName, sepValue, headerValue)}
\arguments{  
  \item{fileName}{name of the CSV file containing te interactions}
  \item{sepValue}{the same as "sep" in read.csv function,it is the value of the field separator character.}
  \item{headerValue}{whether the CSV file has a header or not, TRUE if the file has a header row, FALSE otherwise}
}
\value{  A list containing the interactions. For each gene/protein, the is an entry in the list with "name" 
containing name of the gen/protein and "interactors" containing the list of genes/proteins interacting with it.
}
\author{Kircicegi Korkmaz}
\examples{
	##-interactionList <- readUndirectedInteractionsFromCsv("Arabidopsis_BioGRID-3.1.72.entrezid.csv", " ", FALSE);
}
\keyword{file}
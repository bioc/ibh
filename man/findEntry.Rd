\name{findEntry}
\Rdversion{1.1}
\alias{findEntry}
\title{
Find the index of an entry in an interaction list
}
\description{
Find the index of an entry in an interaction list, used mostly for internal purposes
}
\usage{
findEntry(interactionList, name)
}
\arguments{
  \item{interactionList}{
list of interactions
}
  \item{name}{
name to be searched in the interaction list
}
}
\value{
an integer that is the index of name in the interaction list
}
\examples{
	require(simpIntLists)
	data(ArabidopsisBioGRIDInteractionEntrezId)
	findEntry(ArabidopsisBioGRIDInteractionEntrezId, 832179)
  }
\keyword{ file }


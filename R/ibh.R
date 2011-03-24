
findEntry <- function(interactionList, name){
  which(sapply(interactionList, "[[", "name")  %in% name)
}

ibh <- function(interactionList , geneList){
  geneListAsVector <- unlist(geneList)
  entries <- findEntry( interactionList, geneListAsVector)
  tmp <- unlist(sapply(interactionList[entries], "[[", "interactors"))
  sum(tmp %in% geneListAsVector) / (length(geneList)^2);
}

ibhBioGRID <- function(geneList, organism, idType='EntrezId'){
	interactionList <- findInteractionList(organism, idType);
	ibh(interactionList, geneList);
}

ibhClusterEval <- function(cluster, allGenesList, interactionList){
  sapply(1:max(cluster), function(allGenesList,i){ibh(interactionList,allGenesList[cluster==i])},allGenesList)	
}

ibhClusterEvalBioGRID <- function(cluster, allGenesList, organism, idType='EntrezId'){
	sapply(1:max(cluster), function(i, allGenesList,organism, idType)
        {ibhBioGRID(allGenesList[cluster==i], organism, idType)}, allGenesList, organism, idType)	
}

ibhForMultipleGeneLists <- function(interactionList, listofGeneList) {
  sapply(1:length(listofGeneList),function(i,interactionList, listofGeneList) 
                 ibh(interactionList, listofGeneList[[i]]),interactionList, listofGeneList)
}

ibhForMultipleGeneListsBioGRID <- function(listofGeneList, organism, idType='EntrezId') {
	interactionList <- findInteractionList(organism, idType);
	sapply(1:length(listofGeneList),function(i,interactionList, listofGeneList) 
      ibh(interactionList, listofGeneList[[i]]),interactionList, listofGeneList)
}

readDirectedInteractionsFromCsv <- function(fileName, sepValue, headerValue){
  intFromCsv <- read.csv(fileName, sep = sepValue, header=headerValue);
  l <- list(name=as.vector(intFromCsv[1,1]), interactors=as.vector(intFromCsv[1,2]));
  interactionList <- list();
  interactionList[[1]] <- l;
  for (i in 2:length(intFromCsv[,1])) {
  	if ((intFromCsv[i,1] != "-") & (intFromCsv[i,2] != "-")) {
      entryNumber = findEntry(interactionList, as.vector(intFromCsv[i,1]));
      if (entryNumber == 0) {
        interactionList[[length(interactionList) +1]] <- list(name=as.vector(intFromCsv[i,1]), interactors=as.vector(intFromCsv[i,2]));
      }
      else {
        if (!(intFromCsv[i,2] %in% interactionList[[entryNumber]]$interactors)) {
          interactionList[[entryNumber]]$interactors = c(as.vector(interactionList[[entryNumber]]$interactors), as.vector(intFromCsv[i,2]));
   			}
  		}
  	}
  }
  interactionList;
}

readUndirectedInteractionsFromCsv <- function(fileName, sepValue, headerValue){
  intFromCsv <- read.csv(fileName, sep = sepValue, header=headerValue);
  l <- list(name=as.vector(intFromCsv[1,1]), interactors=as.vector(intFromCsv[1,2]));
  interactionList <- list();
  interactionList[[1]] <- l;
  l2 <- list(name=as.vector(intFromCsv[1,2]), interactors=as.vector(intFromCsv[1,1]));
  interactionList[[2]] <- l2;
  for (i in 2:length(intFromCsv[,1])) {
    if ((intFromCsv[i,1] != "-") & (intFromCsv[i,2] != "-")) {
      entryNumber = findEntry(interactionList, as.vector(intFromCsv[i,1]));
      if (entryNumber == 0) {
        interactionList[[length(interactionList) +1]] <- list(name=as.vector(intFromCsv[i,1]), interactors=as.vector(intFromCsv[i,2]));
      }
      else {
        if (!(intFromCsv[i,2] %in% interactionList[[entryNumber]]$interactors)) {
             interactionList[[entryNumber]]$interactors = 
             c(as.vector(interactionList[[entryNumber]]$interactors), as.vector(intFromCsv[i,2]));
        }
      }

      entryNumber = findEntry(interactionList, as.vector(intFromCsv[i,2]));
      if (entryNumber == 0) {
        interactionList[[length(interactionList) +1]] <- list(name=as.vector(intFromCsv[i,2]), interactors=as.vector(intFromCsv[i,1]));
      }
      else {
        if (!(intFromCsv[i,1] %in% interactionList[[entryNumber]]$interactors)) {
           interactionList[[entryNumber]]$interactors = 
           c(as.vector(interactionList[[entryNumber]]$interactors), as.vector(intFromCsv[i,1]));
        }
      }
    }
  }
interactionList;
}


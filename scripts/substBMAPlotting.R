require(ggplot2)
require(reshape2)

modelNames <- c("JC69", "K80", "F81", "HKY85", "TN93", "GTR")
nModels <- length(modelNames)

readAsList <- function(filename) {
    fh <- file(filename)
    thelist <- lapply(strsplit(readLines(fh), "\t")[-1], as.numeric)
    close(fh)

    thelist <- lapply(thelist, function(x) {x[-1]})

    return(thelist)
}

readTableWithoutFirstColumn <- function(filename, ...) {
    df <- read.table(filename, ...)
    df <- df[,-1]
}

loadData <- function(prefix, burninFrac=0.1) {

    # Read in data
    modelList <- readAsList(paste(prefix,"_modelList_1.log", sep=""))
    substListPrint <-  readAsList(paste(prefix,"_subst.list.print_1.log", sep=""))
    substListPointers <- readTableWithoutFirstColumn(paste(prefix,"_subst.pointers.print_1.log", sep=""), header=T)

    # Remove burnin
    modelList <- modelList[-(1:ceiling(burninFrac*length(modelList)))]
    substListPrint <- substListPrint[-(1:ceiling(burninFrac*length(substListPrint)))]
    substListPointers <- substListPointers[-(1:ceiling(burninFrac*dim(substListPointers)[1])),]

    # Assemble into list
    res <- list()
    res$modelList <- modelList
    res$substListPrint <- substListPrint
    res$substListPointers <- substListPointers

    return(res)
}

getModelPosterior <- function(data) {
    nSites <- dim(data$substListPointers)[2]
    nSamples <- dim(data$substListPointers)[1]

    modelFreqs <- matrix(data=0, nrow=nModels, ncol=nSites)
    rownames(modelFreqs) <- modelNames

    for (i in 1:nSamples) {
        for (s in 1:nSites) {
            model <- data$modelList[[i]][which(data$substListPrint[[i]]==data$substListPointers[i,s])] + 1
            modelFreqs[model,s] <- modelFreqs[model,s] + 1
        }
    }

    modelProbs <- apply(modelFreqs, 2, function(col) {return(col/sum(col))})

    return(modelProbs)
}

plotModelPosterior <- function (prefix, main="Per-site substitution model probability") {
    data <- loadData(prefix)
    modelProbs <- getModelPosterior(data)

    df <- melt(modelProbs)
    names(df) <- c("Model", "Site", "Prob")

    p <- ggplot(df) + geom_tile(aes(x=Site, y=Model, fill=Prob)) + ggtitle(main)

    print(p)
}

plotModelPosterior("mammal_468_sdpm1")

#data <- loadData("mammal_468_sdpm1")

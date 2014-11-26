# This command will look in the XXX folder for files ending in `*subreads.fasta`
# and will count the number of times each sequence header occurs. That number is
# used as the coverage for the read

makeCoverageFiles <- function(folder){
  suffix <- ".subreads.fasta"

  file <- paste(folder, "/", folder, suffix, sep="")
  data <- scan(file, what="", quiet=T)

  headers <- data[grepl("^>", data)]
  headers <- substr(headers, 2, length(headers))

  id <- gsub(".*/(\\d*)/.*", "\\1", headers)
  start <- as.numeric(gsub(".*/(\\d*)_\\d*", "\\1", headers))
  end <- as.numeric(gsub(".*/\\d*_(\\d*)", "\\1", headers))

  freq <- table(id)
  length <- aggregate(end-start, by=list(id), sum)

  data <- matrix(rep(NA, 2*length(freq)), ncol=2)
  colnames(data) <- c("freq", "length")
  rownames(data) <- names(freq)
  data[,"freq"] <- freq[]
  data[,"length"] <- length$x
  write.table(data[order(as.numeric(rownames(data))),], file=paste(folder, "/", folder, ".coverage", sep=""), quote=F)
    return(folder)
}


# This command will look in the XXX folder to see whether the appropriate
# `*.coverage` files are present. If they aren't then it will run
# `makeCoverageFiles`. Once those files are found it then reads the names of the
# sequences that were assembled to make a consensus sequence and extracts the
# coverage data for those sequences.

getCoverage <- function(region){
  write(region, "")


  coverage <- read.table(file=paste("subreads.fasta/", region, "/", region, ".coverage", sep=""), header=T)

  file <- paste(region, "/", list.files(region, pattern="v\\d*.mock?.fasta"), sep="")
  fasta <- scan(file, what="", quiet=T)
  seqNames <- fasta[grepl(">", fasta)]
  seqNames <- gsub(".*/(\\d*)/.*", "\\1", seqNames)

  mock.coverage <- coverage[(seqNames),]

  write.table(mock.coverage, file=paste(region, "/", region, ".ccs.coverage", sep=""))
}

# This should be within the Rmd file...
lapply(dir("./", pattern="v\\d*"), getCoverage)






# Let's also calculate the average quality score for all of the trimmed
# sequences as well as the minimum average quality score over a 50 bp window.
# These two functions assume that the quality scores are coming in as a string
# with spaces between the numbers. Also, any number over 72 should be a 72.

getAverageScore <- function(scoreString){
  scoreVector <- as.numeric(unlist(strsplit(scoreString, " ")))
    scoreVector[scoreVector > 72] <- 72
  return(mean(scoreVector))
}

getMinRollingAverage <- function(scoreString, windowSize){
  scoreVector <- as.numeric(unlist(strsplit(scoreString, " ")))
  scoreVector[scoreVector > 72] <- 72
  min.avg <- 0

  if(length(scoreVector) > windowSize){
    avg <- filter(scoreVector, rep(1/windowSize, windowSize), side=1)
    min.avg <- min(avg, na.rm=T)
  }
  return(min.avg)
}


# Here we repor the average and minimum rolling average quality score across
# each sequencing read. The input is folder name and it looks to see whether
# *.mock.qreport is present and/or newer than the *.mock.qual file.

reportAverageScores <- function(folder, windowSize = 50){
  write(folder, "")
  qual <- scan(file=paste(folder, "/", folder, ".mock.qual", sep=""), what="", sep="\n", quiet=T)
  seq.names <- qual[1:length(qual) %% 2 == 1]
  seq.names <- gsub(".*/(\\d*)/.*", "\\1", seq.names)
  seq.scores <- qual[1:length(qual) %% 2 == 0]

  ave.scores <- unlist(lapply(seq.scores, getAverageScore))
  min.scores <- unlist(lapply(seq.scores, getMinRollingAverage, windowSize))

  report <- cbind(ave.scores, min.scores)
  rownames(report) <- seq.names
  colnames(report) <- c("ave", "min")

  write.table(report, paste(folder, "/", folder, ".mock.qreport", sep=""), quote=F)
}

# This should be within the Rmd file...
lapply(dir("./", pattern="v\\d*"), reportAverageScores)





# We need to calculate the mode to find the best start and end positions in the
# reference alignment

getMode <- function(x){
  return(as.numeric(names(sort(table(x), decreasing=T)[1])))
}


# Here we take in a folder/region and make a composite data table that has all
# of the information about alignment, quality, errors, etc

generateComposite <- function(folder){
  write(folder, "")

  #read everything in
  coverage <- read.table(file=paste(folder, "/", folder, ".ccs.coverage", sep=""), header=T, row.names=1)
  mismatches <- read.table(file=paste(folder, "/", folder, ".mismatches", sep=""), header=F, row.names=1)
  aveq <- read.table(file=paste(folder, "/", folder, ".mock.qreport", sep=""), header=F, skip=1)
  error <- read.table(file=paste(folder, "/", folder, ".mock.filter.error.summary", sep=""), header=T, row.names=1)
  summary <- read.table(file=paste(folder, "/", folder, ".mock.filter.summary", sep=""), header=T, row.names=1)

  #remove chimeras
  non.chimeras <- error$numparents==1

  coverage <- coverage[non.chimeras,]
  mismatches <- mismatches[non.chimeras,]
  aveq <- aveq[non.chimeras,]
  error <- error[non.chimeras,]
  summary <- summary[non.chimeras,]

  #fix some column names
  colnames(mismatches) <- c("barcode", "primer")
  aveq <- aveq[,-1]
  colnames(aveq) <- c("aveQ", "minQ")

  mode.start <- getMode(summary$start)
  mode.end <- getMode(summary$end)

  good.start <- summary$start == mode.start  #find sequences that start at correct location in alignment
  good.end <- summary$end == mode.end        #find sequences that end at correct location in alignment
  good.homop <- summary$polymer <= 8         #find sequences with less than or equal to 8 nt
  good.ambig <- summary$ambig == 0           #find sequences with no ambiguous base calls

  reason <- rep("x", n.seqs)
  reason <- if(good.start, , paste0(reason, "s"))
  reason <- if(good.end, , paste0(reason, "e"))
  reason <- if(good.homop, , paste0(reason, "h"))
  reason <- if(good.ambig, , paste0(reason, "n"))
  reason <- gsub("x", "", reason)
  reason[reason == ""] <- "g"

#	good.sequences <- good.start & good.end & good.homop & good.ambig

  #create composite data frame of good sequences
  composite <- cbind(error[good.sequences,], coverage[good.sequences,],
    mismatches[good.sequences,], aveq[good.sequences,],
    summary[good.sequences,], reason)
  write.table(composite, paste(folder, "/", folder, ".composite", sep=""), quote=F, sep="\t")
}

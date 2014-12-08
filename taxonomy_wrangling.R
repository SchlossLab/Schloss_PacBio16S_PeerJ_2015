regions <- c("v4", "v13", "v35", "v15", "v16", "v19")


# this will go through each vector of taxnomic levels and extract
# bootstrap value. It returns whether a taxonomic level had a bootstrap
# valu over the cutoff (80%)
countGoodBootstraps <- function(nameConf, cutoff=80){
	bootstrap <- gsub(".*\\((\\d*)\\).*", "\\1", nameConf)
	bootstrap[bootstrap == "unclassified"] = 0
	return(sum(as.numeric(bootstrap) >= cutoff))
}


# This function will read in a taxonomy file and parse it for the
# sequence name and the depth of the classification that has a
# boostrap value over the cutoff (80%)
getDepths <- function(taxFileName, cutoff=80){
	tax <- scan(taxFileName, what="", quiet =T)
	lines <- 1:length(tax)
	seqNames <- tax[lines %% 2 == 1]
	taxString <- tax[lines %% 2 == 0]
	taxList <- strsplit(taxString, ";")	#take the taxnomy string and split it by the ; into a list 
	depths <- unlist(lapply(taxList, countGoodBootstraps))
	names(depths) <- seqNames
	return(depths)
}

# go thorugh each region and read in the data from each of the taxonomy files and get the
# classification depth for each sequence read
for(r in regions){
	
	# here we're looking at the rdp, gg, and silva classifications for the
	# "observed" data that used UCHIME to call chimeras
	rdpFileName <- paste("analysis/", r, "/", r, ".trim.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy", sep="")
	rdp <- getDepths(rdpFileName)
	
	ggFileName <- paste("analysis/", r, "/", r, ".trim.unique.good.filter.unique.precluster.pick.gg.wang.taxonomy", sep="")
	gg <- getDepths(ggFileName)
	
	silvaFileName <- paste("analysis/", r, "/", r, ".trim.unique.good.filter.unique.precluster.pick.bacteria.wang.taxonomy", sep="")
	silva <- getDepths(silvaFileName)
	
	
	# output to a file ending in *.tax.compare a list of all sequences and the
	# classification depth for each sequence by the database that was used
	write.table(cbind("rdp" = rdp, "gg" = gg, "silva"=silva), file=paste("analysis/", r, "/", r, ".tax.compare", sep=""), quote=F)
	
	
# i'm not sure what the following lines are for...
#	rdpMock <- paste("analysis/", r, "/", r, ".mock.precluster.pds.wang.taxonomy", sep="")
#	rdp <- getDepths(rdpMock)
#	
#	ggMock <- paste("analysis/", r, "/", r, ".mock.precluster.gg.wang.taxonomy", sep="")
#	gg <- getDepths(ggMock)
#	
#	silvaMock <- paste("analysis/", r, "/", r, ".mock.precluster.bacteria.wang.taxonomy", sep="")
#	silva <- getDepths(silvaMock)
#	
#	write.table(cbind("rdp" = rdp, "gg" = gg, "silva"= silva), file=paste("analysis/", r, "/mock.tax.compare", sep=""), quote=F)  
}


# here we want to merge the classification depth that we got from each sequence
# as well as the number of sequences that each sequence represents by library
getDepthByLibrary <- function(region){
	
	# get the frequency data...
	count.file <- paste("analysis/", region, "/", region, ".trim.unique.good.filter.unique.precluster.pick.count_table", sep="")
	count.table <- read.table(file=count.file, header=T, row.names=1)
	
	# get the read depth data...
	depth.file <- paste("analysis/", region, "/", region, ".tax.compare", sep="")
	depth.table <- read.table(file=depth.file, header=T, row.names=1)
	depth.table$rdp <- factor(depth.table$rdp, levels=0:7)
	depth.table$gg <- factor(depth.table$gg, levels=0:7)
	depth.table$silva <- factor(depth.table$silva, levels=0:7)
	
	mock.depth <- depth.table[count.table$mock > 0,]
	human.depth <- depth.table[count.table$human > 0,]
	mouse.depth <- depth.table[count.table$mouse > 0,]
	soil.depth <- depth.table[count.table$soil > 0,]
	
	# calculate the percentages for each database
	mock.gg <- 100*summary(mock.depth$gg)/nrow(mock.depth)
	human.gg <- 100*summary(human.depth$gg)/nrow(human.depth)
	mouse.gg <- 100*summary(mouse.depth$gg)/nrow(mouse.depth)
	soil.gg <- 100*summary(soil.depth$gg)/nrow(soil.depth)
	
	mock.rdp <- 100*summary(mock.depth$rdp)/nrow(mock.depth)
	human.rdp <- 100*summary(human.depth$rdp)/nrow(human.depth)
	mouse.rdp <- 100*summary(mouse.depth$rdp)/nrow(mouse.depth)
	soil.rdp <- 100*summary(soil.depth$rdp)/nrow(soil.depth)
	
	mock.silva <- 100*summary(mock.depth$silva)/nrow(mock.depth)
	human.silva <- 100*summary(human.depth$silva)/nrow(human.depth)
	mouse.silva <- 100*summary(mouse.depth$silva)/nrow(mouse.depth)
	soil.silva <- 100*summary(soil.depth$silva)/nrow(soil.depth)
	
	return(rbind(mock.rdp, human.rdp, mouse.rdp, soil.rdp, mock.gg, human.gg, mouse.gg, soil.gg, mock.silva, human.silva, mouse.silva, soil.silva))
}

# for each region get the percentage of sequences in each library that
# classified to each taxonomic level for each databaes considered
composite <- data.frame(matrix(rep(0, 8*6*12), ncol=8))
colnames(composite) <- 0:7
composite[1:12,] <- getDepthByLibrary("v4");
composite[13:24,] <- getDepthByLibrary("v35")
composite[25:36,] <- getDepthByLibrary("v13")
composite[37:48,] <- getDepthByLibrary("v15")
composite[49:60,] <- getDepthByLibrary("v16")
composite[61:72,] <- getDepthByLibrary("v19")

# format the final output file
composite$region <- c(rep("v4", 12), rep("v35", 12), rep("v13", 12), rep("v15", 12), rep("v16", 12), rep("v19", 12))
composite$database <- rep(c(rep("rdp", 4), rep("gg", 4), rep("silva", 4)), 6)
composite$sample <- rep(c("mock", "human", "mouse", "soil"), 18)
composite$total <- composite[,"6"] + composite[,"7"] #make a total column that has the % of genus and species-level names
write.table(file="taxonomy.depth.analysis", composite, quote=F)

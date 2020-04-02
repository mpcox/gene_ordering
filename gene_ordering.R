# simulate random deletion of genes in the yeast genome
# and identify shared gene ordering

# variables
iterations       <- 100
n.chr            <- 8
n.anc.genes      <- 4943
n.postdupl.genes <- 5500
filename         <- "outfile.txt"

# make header line
cat("n.shared\tshared.gene.pairs\tsim1.duplicates\tsim2.duplicates\tsim1.gen1.pairs\tsim1.gen2.pairs\tsim2.gen1.pairs\tsim2.gen2.pairs\n", file="outfile.txt", append=T)

# make genome and duplicate
genome <- seq(1, n.anc.genes - n.chr)
genome <- c(genome, genome)

# run iterations
shared.counts <- vector()
for( x in 1:iterations ){

	# generate simulated dataset 1
	# randomly remove genes until n.postdupl.genes remain
	sim1 <- genome
	sim1[sample(1:length(genome), length(genome)-n.postdupl.genes)] <- NA
	
	# break into chromosomes
	sim1.gen1 <- sim1[1:n.anc.genes]
	sim1.gen2 <- sim1[(n.anc.genes+1):(2*n.anc.genes)]
	
	# identify duplicate genes
	sim1.duplicates <- na.omit(intersect(sim1.gen1, sim1.gen2))
	
	# simplify by removing NA states
	sim1.gen1 <- sim1.gen1[!is.na(sim1.gen1)]
	sim1.gen2 <- sim1.gen2[!is.na(sim1.gen2)]
	
	# make lists of adjacent gene pairs
	sim1.gen1.pairs <- vector()
	sim1.gen2.pairs <- vector()
	for( i in 1:(length(sim1.gen1)-1) ){
		sim1.gen1.pairs <- c(sim1.gen1.pairs, paste0(sim1.gen1[i], "_", sim1.gen1[i+1]) )
	}
	for( i in 1:(length(sim1.gen2)-1) ){
		sim1.gen2.pairs <- c(sim1.gen2.pairs, paste0(sim1.gen2[i], "_", sim1.gen2[i+1]) )
	}
	sim1.pairs <- c(sim1.gen1.pairs, sim1.gen2.pairs)
	
	# remove any pair containing a duplicate gene
	list <- strsplit(sim1.pairs, "_")
	df <- data.frame(matrix(unlist(list), nrow=length(list), byrow=T), stringsAsFactors=F)
	sim1.pairs[(which( df[,1] %in% sim1.duplicates ))] <- NA
	sim1.pairs[(which( df[,2] %in% sim1.duplicates ))] <- NA
	
	# remove new NA states
	sim1.pairs <- sim1.pairs[!is.na(sim1.pairs)]
	
	
	# generate simulated dataset 2
	# randomly remove genes until n.postdupl.genes remain
	sim2 <- genome
	sim2[sample(1:length(genome), length(genome)-n.postdupl.genes)] <- NA
	
	# break into chromosomes
	sim2.gen1 <- sim2[1:n.anc.genes]
	sim2.gen2 <- sim2[(n.anc.genes+1):(2*n.anc.genes)]
	
	# identify duplicate genes
	sim2.duplicates <- na.omit(intersect(sim2.gen1, sim2.gen2))
	
	# simplify by removing NA states
	sim2.gen1 <- sim2.gen1[!is.na(sim2.gen1)]
	sim2.gen2 <- sim2.gen2[!is.na(sim2.gen2)]
	
	# make lists of adjacent gene pairs
	sim2.gen1.pairs <- vector()
	sim2.gen2.pairs <- vector()
	for( i in 1:(length(sim2.gen1)-1) ){
		sim2.gen1.pairs <- c(sim2.gen1.pairs, paste0(sim2.gen1[i], "_", sim2.gen1[i+1]) )
	}
	for( i in 1:(length(sim2.gen2)-1) ){
		sim2.gen2.pairs <- c(sim2.gen2.pairs, paste0(sim2.gen2[i], "_", sim2.gen2[i+1]) )
	}
	sim2.pairs <- c(sim2.gen1.pairs, sim2.gen2.pairs)
	
	# remove any pair containing a duplicate gene
	list <- strsplit(sim2.pairs, "_")
	df <- data.frame(matrix(unlist(list), nrow=length(list), byrow=T), stringsAsFactors=F)
	sim2.pairs[(which( df[,1] %in% sim2.duplicates ))] <- NA
	sim2.pairs[(which( df[,2] %in% sim2.duplicates ))] <- NA
	
	# remove new NA states
	sim2.pairs <- sim2.pairs[!is.na(sim2.pairs)]
	
	
	# find shared pairs
	shared <- intersect(sim1.pairs,sim2.pairs)
	n.shared <- length(shared)
	
	# write information to file
	if( length(shared) == 0 ){
		shared <- "NA"
	}
	cat(c(n.shared, "\t", shared, "\t"), file=filename, append=T)
	cat(c(sim1.duplicates, "\t", sim2.duplicates, "\t"), file=filename, append=T)
	cat(c(sim1.gen1.pairs, "\t", sim1.gen2.pairs, "\t", sim2.gen1.pairs, "\t", sim2.gen2.pairs, "\n"), file=filename, append=T)
	shared.counts <- c(shared.counts, n.shared)

}

hist(shared.counts)

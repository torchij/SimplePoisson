args = commandArgs(trailingOnly=TRUE)

if (length(args) != 3){
    cat("Usage: get_ppois_from_bam.r bam_file bed_file cpu_number\n")
    quit(save="no")
} else {
    library(rtracklayer)
    library(data.table)
    #library(seqbias)
    #library(cn.mops)
    #library(vcd)
    library(GenomicAlignments)
}

#bam_file <- args[1] # test.bam
#bed_file <- args[2] # test.bed
#con_file <- args[3] # cont.bed
#cpu_numb <- args[4] # cpu numb

# inputs test 3
#setwd("/Users/jtorchia/git/torchij/SimplePoisson/")
#bam_file <- "hcc1187_filtered.sorted.regions_v.b.2.1.bam"     # input bam file (**NB** PREFILTERED TRANS-ONLY BAM)
#bed_file <- "regions_v.b.2.1.bed"
#con_file <- "cont.bed"                              # control probe bed file
cpu_numb <- 2                                       # parallelize the huge bam counts step

# sampling inputs
num_reps <- 1000

# set output name
samp <- gsub("\\.bam", "", basename(bam_file))

poissonness_test <- function(counts) {

    # based on Hoaglin, 1980: “A poissonness plot”
    n=length(counts)
    x=table(counts)
    k=as.numeric(names(x))
    f=c(log(x) + lfactorial(k))
    r=n*(1-cor(k,log(x)+lfactorial(k))^2)
    ppois_test <- list(k, f, r)
    return(ppois_test)

}

##################################################
##### Step 1: Calculate signal from bam file  ####
##################################################

# read in bam file
galp <- readGAlignmentPairs(bam_file, use.names=TRUE, strandMode=1)

# get trans reads
galp <- galp[seqnames(first(galp)) != seqnames(second(galp))]

# read in bed file with probes
bed <- read.table(bed_file, header=FALSE)
colnames(bed) <- c("chrom", "start", "end", "tlocID", "score", "strand") # for region file

# aggregate concatenate to get chromosome pairs
pre <- unique(bed[c(1:4)])
key <- aggregate(cbind(chrom, start, end) ~ tlocID, data=unique(bed[c(1:4)]), FUN = c)

# for each row, do stuff
padding=10000
result <- data.frame(matrix(ncol=8, nrow=0))
resNam <- c("tlocID", "chrA", "startA", "endA", "chrB", "startB", "endB", "support")
colnames(result) <- resNam
for (row in 1:nrow(key)) {

    tloc <- key[row, "tlocID"]
    chrA <- key[row, "chrom"][1]
    chrB <- key[row, "chrom"][2]
    staA <- as.numeric(key[row, "start"][1]) - padding
    staB <- as.numeric(key[row, "start"][2]) - padding
    endA <- as.numeric(key[row, "end"][1]) + padding
    endB <- as.numeric(key[row, "end"][2]) + padding
    a <- length(galp[ (seqnames(first(galp)) == chrA & (start(first(galp)) > staA & start(first(galp)) < endA)) & seqnames(second(galp)) == chrB ])
    b <- length(galp[ (seqnames(first(galp)) == chrB & (start(first(galp)) > staB & start(first(galp)) < endB)) & seqnames(second(galp)) == chrA ])
    c <- length(galp[ (seqnames(second(galp)) == chrA & (start(second(galp)) > staA & start(second(galp)) < endA)) & seqnames(first(galp)) == chrB ])
    d <- length(galp[ (seqnames(second(galp)) == chrB & (start(second(galp)) > staB & start(second(galp)) < endB)) & seqnames(first(galp)) == chrA ])
    support <- a + b + c + d
    #support <- length(galp[ (seqnames(first(galp)) == chrA & seqnames(second(galp)) == chrB) | (seqnames(first(galp)) == chrB & seqnames(second(galp)) == chrA) ])
    newRow <- data.frame(tlocID=tloc, chrA=chrA, startA=staA, endA=endA, chrB=chrB, startB=staB, endB=endB, support=support)
    result <- rbind(result, newRow)

}

##################################################
##### Step 2: Estimate lambda from bam file  #####
##################################################

# we use bedtools to get the counts per probe interval
# side note: there may be an r-native way to do this,
# but bedtools does it so efficiently...
cmd <- paste0("bedtools multicov -bams ", bam_file, " -bed ", con_file)

# fread has a nice method to take output of system command
df_con <- fread(cmd=cmd)

# rename columns
colnames(df_con) <- c("chrom", "start", "end", "probe", "score", "strand", "counts")

# sample n rows from control (17 is the number of probes per tloc)
counts_sum <- replicate(num_reps, sum(df_con[sample(nrow(df_con), 34), ]$counts))
lambda <- median(counts_sum) + 1

# make a table for output
lambda_df <- data.frame(lambda=lambda, counts_avg=mean(counts_sum), counts_std=sd(counts_sum))

######################################################################################
##### Step 3: Test poissonness of background and plot the empirical distribution  ####
######################################################################################

# poissonness test for background and foreground
bg <- poissonness_test(counts_sum)
k=bg[[1]]; f=bg[[2]]; r=bg[[3]]
fg <- poissonness_test(counts_sum)
k2=fg[[1]]; f2=fg[[2]]; r2=fg[[3]]
png(paste0(samp, ".gf.png"))
    plot(k, f, col="black", main=paste0("poissonness test (bg) - ", as.character(signif(r, 4))))
    points(k2, f2, col="red")
dev.off()

# plot empirical distribution of background
png(paste0(samp, ".bg.hist.png"))
    hist(counts_sum, breaks=50, main="background counts distribution")
dev.off()

#########################################################
##### Step 4: Apply poisson model to signal, lambda #####
#########################################################

# calculate ppois (lambda is the average counts aggregated over a control sample)
result$pvalue <- ppois(q=result$support, lambda=lambda, lower.tail=FALSE, log=FALSE)

# write to table
write.table(result, file=paste0(samp, ".results.txt"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(lambda_df, file=paste0(samp, ".lambda.txt"), sep="\t", row.names=FALSE, quote=FALSE)
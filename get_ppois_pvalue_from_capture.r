args = commandArgs(trailingOnly=TRUE)

if (length(args) != 3){
    cat("Usage: get_ppois_from_bam.r bam_file bed_file cpu_number\n")
    quit(save="no")
} else {
    library(rtracklayer)
    library(data.table)
    library(seqbias)
    library(cn.mops)
    library(vcd)
}

#bam_file <- args[1] # test.bam
#bed_file <- args[2] # test.bed
#con_file <- args[3] # cont.bed
#cpu_numb <- args[4] # cpu numb

# inputs test 1
#setwd("/Users/jtorchia/git/torchij/SimplePoisson/")
#bam_file <- "test.bam"      # input bam file (**NB** PREFILTERED TRANS-ONLY BAM)
#bed_file <- "test.bed"      # target probe bed file
#con_file <- "cont.bed"      # control probe bed file
#cpu_numb <- 2               # parallelize the huge bam counts step

# inputs test 2
setwd("/Users/jtorchia/git/torchij/SimplePoisson")
bam_file <- "test.bam"                              # input bam file
bed_file <- "test.bed"                              # target probe bed file
con_file <- "cont.bed"                              # control probe bed file
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
##### Step 1: Estimate lambda from bam file  #####
##### (with control probes)                  #####
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
lambda <- median(counts_sum)

# make a table for output
lambda_df <- data.frame(lambda=lambda, counts_avg=mean(counts_sum), counts_std=sd(counts_sum))

##################################################
##### Step 2: Calculate normalized signal    #####
##################################################

# get counts per interval for the signal probes
cmd <- paste0("bedtools multicov -bams ", bam_file, " -bed ", bed_file)

# fread has a nice method to take output of system command
df <- fread(cmd=cmd)

# rename columns
colnames(df) <- c("chrom", "start", "end", "size", "strand", "tloc", "probe", "counts")

# aggregate by tloc
df_agg <- aggregate(counts ~ tloc, data = df, FUN = sum, na.rm = TRUE)

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
df_agg$pvalue <- ppois(q=df_agg$counts, lambda=lambda, lower.tail=FALSE, log=FALSE)

# write to table
write.table(df_agg, file=paste0(samp, ".results.txt"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(lambda_df, file=paste0(samp, ".lambda.txt"), sep="\t", row.names=FALSE, quote=FALSE)
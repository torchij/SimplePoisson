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

bam_file <- args[1] # test.bam
bed_file <- args[2] # test.bed
cpu_numb <- args[3]

# inputs
setwd("/Users/jtorchia/git/torchij/SimplePoisson/")
bam_file <- "test.bam"      # input bam file
bed_file <- "test.counts"   # target probe bed file with read counts per probe
cpu_numb <- 2             # parallelize the huge bam counts step

# sampling inputs
bin_size <- 1000
num_intv <- 1000
num_reps <- 1000

# set output name
samp <- gsub("\\.bam", "", basename(bam_file))

estimate_lambda <- function(gr, bin_size, num_intv, num_reps) {
   
    int_size <- bin_size * num_intv

    ### Method: sample n rows (num_intv) of binned counts
    ### The goal is to select random regions of the genome,
    ### sum the reads in that region, and define that as the
    ### background rate
    counts_sum <- replicate(num_reps, sum(sample(gr$sample, size=num_intv)))
    counts_avg <- mean(counts_sum)
    counts_std <- sd(counts_sum)
    lambda <- counts_avg / int_size * 1000 # we set as a rate per kb
    lambda_df <- data.frame(lambda=lambda, counts_avg=counts_avg, counts_std=counts_std)
    results <- list(counts_sum, lambda_df)

    # return list of counts stats and the lambda
    return(results)

}

poissonness_test <- function(counts_sum) {

    # based on Hoaglin, 1980: “A poissonness plot”
    n=length(counts_sum)
    x=table(counts_sum)
    k=as.numeric(names(x))
    f=c(log(x) + lfactorial(k))
    r=n*(1-cor(k,log(x)+lfactorial(k))^2)
    ppois_test <- list(k, f, r)
    return(ppois_test)

}

##################################################
##### Step 1: Estimate lambda from bam file  #####
##################################################

# get counts by bin size
gr <- getReadCountsFromBAM(bam_file, sampleNames="sample", WL=bin_size, parallel=cpu_numb)

# Estimate lambda
results <- estimate_lambda(gr, bin_size, num_intv, num_reps)
lambda <- results[[2]]$lambda
counts_sum <- results[[1]]

##################################################
##### Step 2: Calculate normalized signal    #####
##################################################

# we use bedtools to get the counts per probe interval
# side note: there may be an r-native way to do this,
# but bedtools does it so efficiently...
cmd <- paste0("bedtools multicov -bams ", bam_file, " -bed ", bed_file)

# fread has a nice method to take output of system command
df <- fread(cmd=cmd)

# rename columns
colnames(df) <- c("chrom", "start", "end", "tloc", "score", "strand", "counts")

# add size of interval
df$size <- df$end - df$start

# aggregate mean of coverage and control by translocation
df_agg <- aggregate(cbind(counts, size) ~ tloc, data = df, FUN = sum, na.rm = TRUE)

# normalize counts by size
df_agg$counts_kb <- df_agg$counts / df_agg$size * 1000

######################################################################################
##### Step 3: Test poissonness of background and plot the empirical distribution  ####
######################################################################################

# poissonness test for background and foreground
bg <- poissonness_test(counts_sum)
k=bg[[1]]; f=bg[[2]]; r=bg[[3]]
fg <- poissonness_test(df$counts)
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
df_agg$pvalue <- ppois(q=df_agg$counts_kb, lambda=lambda, lower.tail=FALSE, log=FALSE)

# write to table
write.table(df_agg, file=paste0(samp, ".results.txt"), sep="\t", row.names=FALSE)
write.table(results[[2]], file=paste0(samp, ".lambda.txt"), sep="\t", row.names=FALSE)
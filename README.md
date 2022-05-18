# SimplePoisson

SimplePoisson is a basic r script which takes in a bam file from a capture experiment, and a bed file of the capture probe and targeted regions, and outputs a poisson-derived probability that the capture regions are enriched over background. In a way, its kind of like a super simplified peak-calling script, but only when data follows a poisson distribution.

The script does four basic things:

1) Estimate background lambda by sampling read counts from random regions of the bam file

2) Aggregates counts over the targeted regions

3) Outputs some basic plots of the background distribution and a test for poissonness

4) Generates a p-value based on a basic poisson distribution

## Prerequisites

The following R libraries are required:

```bash
    rtracklayer
    data.table
    seqbias
    cn.mops
    vcd
```

## Installation

```bash
git clone https://github.com/torchij/SimplePoisson.git
```

## Usage

```bash
bam=test.bam # the capture bam file
bed=test.bed # the capture probes (chr, start, end, region_name, probe_name)
cpus=4       # available cpus for parallelization

# running
Rscript get_ppois_pvalue_from_capture.r $bam $bed $cpus
```

## License
[MIT](https://choosealicense.com/licenses/mit/)
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(IRanges))


args <- commandArgs(T)
frag <- args[1]

ofname <- file.path(dirname(frag), paste0(gsub("\\..+", "", basename(frag)), ".FSD.csv"))


# read chromosome arms
arms <- "/public/home/kai/projects/fragmentome/fragnorama/sources/hg19_arms.txt"
# arms <- "\\\\172.16.11.242/test_data/fragmentome/sources/hg19_arms.txt"
arms <- fread(arms)


# read fragments
# frag <- "\\\\172.16.11.242/test_data/fragmentome/results/1_benigns/MRD220814-505/frag/MRD220814-505_tumor.step2_filt.bed"
frag <- fread(frag,
              col.names = c("chr", "start", "end", "mapq", "fraggc"))

frag[, "length":=(end-start)]


# annotate chromosome arm into each fragments
setkey(frag, chr, start, end)
setkey(arms, chrom, start, end)

chroms <- paste0("chr",c(1:22, "X"))

frag <- foverlaps(frag[chr %in% chroms], 
                  arms, 
                  type = "within",
                  nomatch = NULL)

frag <- frag[, -c("start", "end")]



# counts fragments in each arm bin
bins <- data.table(bin_start = seq(60, 395, 5),
                   bin_end = seq(65, 400, 5))

concat_bins <- list()

i <- 1
for (chrom in chroms){
    for (iarm in c("p", "q")){
        arm_bin <- IRanges(start = bins$bin_start, 
                           end = bins$bin_end,
                           chrom = chrom,
                           arm = iarm)
        
        arm_frag <- frag[chr == chrom & arm == iarm, ]
        arm_frag <- IRanges(start = arm_frag$length,
                            end = arm_frag$length)
        
        counts <- countOverlaps(arm_bin, arm_frag)
        arm_bin <- as.data.table(arm_bin)
        arm_bin[, "count":=counts]
        concat_bins[[i]] <- arm_bin
        i <- i + 1
    }
}

bins <- do.call(rbind, concat_bins)

fwrite(bins, file = ofname, sep = ",")


















